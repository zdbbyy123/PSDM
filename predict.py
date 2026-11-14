"""Apply the frozen PSDM ensemble to VCF variants and independently prepared annotations."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import zipfile

import numpy as np
import lightgbm as lgb

DEFAULT_MODEL = Path(__file__).with_name('PSDM.model.zip')
PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def text_file(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8-sig', newline='') if path.suffix == '.gz' else path.open(encoding='utf-8-sig', newline='')


def variant_key(chrom, pos, ref, alt):
    chrom = str(chrom).strip()
    chrom = chrom[3:] if chrom.lower().startswith('chr') else chrom
    chrom = chrom.upper()
    if chrom == 'MT':
        chrom = 'M'
    if chrom not in {str(i) for i in range(1, 23)} | {'X', 'Y', 'M'}:
        raise ValueError(f'Unsupported chromosome: {chrom}')
    if not str(pos).isdigit() or int(pos) < 1:
        raise ValueError(f'Invalid 1-based position: {pos}')
    ref, alt = str(ref).upper(), str(alt).upper()
    if ref not in 'ACGT' or alt not in 'ACGT' or len(ref) != 1 or len(alt) != 1 or ref == alt:
        raise ValueError('Each input record must be one biallelic A/C/G/T SNV with REF != ALT.')
    return f'chr{chrom}:{int(pos)}:{ref}>{alt}'


def read_vcf(path):
    records, keys = [], set()
    format_seen = reference_seen = header_seen = False
    with text_file(path) as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('##fileformat=VCFv4.'):
                format_seen = True
            elif line.startswith('##reference='):
                if line.split('=', 1)[1] not in ('GRCh38', 'hg38'):
                    raise ValueError('VCF must explicitly declare ##reference=GRCh38 (or hg38). No liftover is performed.')
                reference_seen = True
            elif line.startswith('#CHROM\t'):
                if header_seen or line.split('\t')[:8] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']:
                    raise ValueError('Invalid or repeated VCF header.')
                header_seen = True
            elif line.startswith('#'):
                continue
            else:
                if not (format_seen and reference_seen and header_seen):
                    raise ValueError('VCF requires fileformat, GRCh38 reference, and #CHROM headers before records.')
                cells = line.split('\t')
                if len(cells) < 8:
                    raise ValueError('VCF record has fewer than eight fields.')
                key = variant_key(cells[0], cells[1], cells[3], cells[4])
                if key in keys:
                    raise ValueError(f'Duplicate VCF variant: {key}')
                keys.add(key)
                records.append({'variant_id': key, 'chrom': key.split(':')[0], 'pos': int(cells[1]),
                                'ref': cells[3].upper(), 'alt': cells[4].upper(), 'input_id': cells[2]})
    if not records:
        raise ValueError('VCF contains no variants.')
    return records


class PSDM:
    """One frozen ensemble: five fold-specific preprocessors and 350 native learners."""

    def __init__(self, model_path=DEFAULT_MODEL):
        self.path = Path(model_path).resolve()
        self._learners = None
        with zipfile.ZipFile(self.path) as archive:
            names = archive.namelist()
            if len(set(names)) != len(names):
                raise ValueError('Model archive has duplicate entries.')
            self.spec = json.loads(archive.read('model.json'))
            if self.spec.get('format') != 'PSDM-inference-1' or self.spec.get('genome_build') != 'GRCh38':
                raise ValueError('Unsupported model format or reference build.')
            self.features = self.spec['features']
            self.names = [item['name'] for item in self.features]
            if len(self.names) != 260 or len(set(self.names)) != 260:
                raise ValueError('The final model requires 260 distinct raw annotations.')
            schema_bytes = json.dumps(self.names, ensure_ascii=True, separators=(',', ':')).encode()
            self.schema_sha256 = hashlib.sha256(schema_bytes).hexdigest()
            if self.schema_sha256 != self.spec['schema_sha256']:
                raise ValueError('Feature schema checksum mismatch.')
            expected = self.spec['files_sha256']
            required = {f'fold_{f}/learner_{b:03d}.txt.gz' for f in range(1, 6) for b in range(70)}
            required |= {f'fold_{f}/preprocessing.npz' for f in range(1, 6)}
            if set(expected) != required or set(names) != required | {'model.json'}:
                raise ValueError('Model must contain exactly 350 learners and five preprocessors.')
            for name, expected_hash in expected.items():
                if hashlib.sha256(archive.read(name)).hexdigest() != expected_hash:
                    raise ValueError(f'Model checksum mismatch: {name}')
            self.preprocessors = []
            for fold in range(1, 6):
                with np.load(io.BytesIO(archive.read(f'fold_{fold}/preprocessing.npz')), allow_pickle=False) as source:
                    prep = dict(source)
                if set(prep) != {'input_indices', 'raw_median', 'mean', 'scale', 'core_median'}:
                    raise ValueError('Unexpected preprocessing parameters.')
                indices = prep['input_indices']
                if indices.shape != (231,) or indices.dtype.kind not in 'iu' or len(set(indices)) != 231 or indices.min() < 0 or indices.max() >= 260:
                    raise ValueError('Invalid fold feature mapping.')
                for name in ('raw_median', 'mean', 'scale'):
                    if prep[name].shape != (96,) or not np.isfinite(prep[name]).all():
                        raise ValueError(f'Invalid preprocessing array: {name}')
                if (prep['scale'] <= 0).any() or prep['core_median'].shape != (231,) or not np.isfinite(prep['core_median']).all():
                    raise ValueError('Invalid scale or final imputation parameters.')
                self.preprocessors.append(prep)

    def _load_learners(self):
        if self._learners is None:
            self._learners = []
            with zipfile.ZipFile(self.path) as archive:
                for fold in range(1, 6):
                    learners = []
                    for bag in range(70):
                        raw = gzip.decompress(archive.read(f'fold_{fold}/learner_{bag:03d}.txt.gz'))
                        # Match gzip.open(..., 'rt'): LightGBM tree offsets require normalized newlines.
                        with io.TextIOWrapper(io.BytesIO(raw), encoding='utf-8', newline=None) as handle:
                            text = handle.read()
                        learner = lgb.Booster(model_str=text)
                        if learner.num_feature() != 231:
                            raise ValueError('Learner feature dimension differs from the frozen model.')
                        learners.append(learner)
                    self._learners.append(learners)

    def predict(self, annotations, feature_names):
        """Predict from a numeric (variants x 260) matrix; column names must match exactly."""
        if list(feature_names) != self.names:
            raise ValueError('Feature names/order differ from model.names; do not infer annotation order.')
        values = np.asarray(annotations, dtype=np.float64)
        if values.ndim != 2 or values.shape[1] != 260 or not len(values) or np.isinf(values).any():
            raise ValueError('Expected a nonempty matrix with 260 numeric columns; NaN is allowed, infinity is not.')
        # A whole uncomputed block is not an acceptable substitute for annotation.
        for block in ('ifeatureomega', 'ftrcool_anf', 'gpn_msa', 'fixed'):
            indices = [i for i, item in enumerate(self.features) if item['block'] == block]
            if np.isnan(values[:, indices]).all(axis=1).any():
                raise ValueError(f'An entire annotation block is missing for at least one variant: {block}')
        self._load_learners()
        fold_predictions = []
        for prep, learners in zip(self.preprocessors, self._learners):
            x = values[:, prep['input_indices']].copy()
            first = x[:, :96]
            x[:, :96] = (np.where(np.isnan(first), prep['raw_median'], first) - prep['mean']) / prep['scale']
            x = np.where(np.isnan(x), prep['core_median'], x)
            if not np.isfinite(x).all():
                raise ValueError('Nonfinite transformed annotations.')
            probabilities = np.column_stack([learner.predict(x, num_threads=1) for learner in learners])
            fold_predictions.append(probabilities)
        # Preserve the original (fold, variant, learner) averaging order.
        scores = np.stack(fold_predictions).mean(axis=(0, 2))
        if not np.isfinite(scores).all() or ((scores < 0) | (scores > 1)).any():
            raise ValueError('Invalid predicted probabilities.')
        return scores

    def read_annotations(self, path, records):
        """Join raw annotation TSV to VCF by chromosome, position, REF and ALT."""
        wanted = {row['variant_id']: i for i, row in enumerate(records)}
        values = np.full((len(records), 260), np.nan)
        seen, found = set(), set()
        build = schema = None
        with text_file(path) as handle:
            header = None
            for line in handle:
                if line.startswith('##reference='):
                    build = line.strip().split('=', 1)[1]
                elif line.startswith('##PSDM_SCHEMA_SHA256='):
                    schema = line.strip().split('=', 1)[1]
                elif line.startswith('##') or not line.strip():
                    continue
                else:
                    header = next(csv.reader([line], delimiter='\t'))
                    break
            if build != 'GRCh38' or schema != self.schema_sha256:
                raise ValueError('Annotations require matching ##reference=GRCh38 and ##PSDM_SCHEMA_SHA256 headers; see --describe.')
            if header is None or len(header) != len(set(header)):
                raise ValueError('Missing or duplicate annotation column names.')
            required = ['chrom', 'pos', 'ref', 'alt'] + self.names
            if not set(required).issubset(header):
                missing = [name for name in required if name not in header]
                raise ValueError('Required annotation columns absent: ' + ', '.join(missing[:8]))
            indices = [header.index(name) for name in required]
            for row in csv.reader(handle, delimiter='\t'):
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError('Annotation row width differs from its header.')
                cells = [row[i] for i in indices]
                key = variant_key(*cells[:4])
                if key in seen:
                    raise ValueError(f'Duplicate annotation variant: {key}')
                seen.add(key)
                if key not in wanted:
                    continue
                parsed = []
                for value in cells[4:]:
                    value = value.strip()
                    number = float('nan') if value.lower() in ('', '.', 'na', 'nan') else float(value)
                    if math.isinf(number):
                        raise ValueError(f'Infinite annotation for {key}.')
                    parsed.append(number)
                values[wanted[key]] = parsed
                found.add(key)
        missing = set(wanted) - found
        if missing:
            raise ValueError(f'{len(missing)} VCF variant(s) lack an annotation row; first: {sorted(missing)[0]}')
        return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, default=DEFAULT_MODEL)
    parser.add_argument('--vcf', type=Path, help='GRCh38 biallelic SNV VCF, optionally gzipped')
    parser.add_argument('--annotations', type=Path, help='Externally prepared, schema-matched raw annotation TSV or TSV.gz')
    parser.add_argument('--output', type=Path, help='New prediction TSV; existing files are not overwritten')
    parser.add_argument('--describe', action='store_true', help='Print the annotation contract and feature schema as JSON')
    args = parser.parse_args()
    if not args.describe and (args.vcf is None or args.annotations is None or args.output is None):
        parser.error('--vcf, --annotations and --output are required. Raw VCF alone cannot supply PSDM annotations.')
    try:
        model = PSDM(args.model)
        if args.describe:
            print(json.dumps({key: model.spec[key] for key in ('format', 'genome_build', 'schema_sha256', 'annotation_contract', 'features')}, indent=2))
            return
        output = args.output.expanduser().resolve()
        if output.exists():
            raise FileExistsError(f'Output already exists: {output}')
        for protected in (PACKAGE_ROOT/'data', PACKAGE_ROOT/'code'):
            if output == protected or protected in output.parents:
                raise ValueError('Choose an output path outside data/ and code/.')
        if output in {args.vcf.resolve(), args.annotations.resolve(), args.model.resolve()}:
            raise ValueError('Output cannot overwrite an input or model.')
        records = read_vcf(args.vcf)
        values = model.read_annotations(args.annotations, records)
        scores = model.predict(values, model.names)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['variant_id', 'chrom', 'pos', 'ref', 'alt', 'input_id', 'psdm_score', 'predicted_label'])
            for row, score in zip(records, scores):
                writer.writerow([row[key] for key in ('variant_id', 'chrom', 'pos', 'ref', 'alt', 'input_id')] + [format(float(score), '.17g'), int(score >= .5)])
        print(f'Predicted {len(records)} variants using all 350 frozen learners: {output}')
    except (ValueError, OSError, KeyError, zipfile.BadZipFile, json.JSONDecodeError, lgb.basic.LightGBMError) as error:
        parser.exit(2, f'Error: {error}\n')


if __name__ == '__main__':
    main()
