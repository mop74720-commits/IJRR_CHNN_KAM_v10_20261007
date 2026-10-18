"""Hash the local scientific closure without importing/re-running experiments."""
from pathlib import Path
import hashlib, importlib.metadata, json, platform

HERE = Path(__file__).resolve().parent


def write(name, value):
    text = json.dumps(value, indent=2, ensure_ascii=False) + '\n'
    (HERE/name).write_bytes(text.replace('\n', '\r\n').encode('utf-8'))


def main():
    versions = {'python':platform.python_version()}
    for name in ['numpy','scipy','numba','matplotlib','threadpoolctl','PyMuPDF']:
        versions[name] = importlib.metadata.version(name)
    protocol = json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    for name, expected in protocol['input_hashes'].items():
        actual = hashlib.sha256((HERE/'inputs'/name).read_bytes()).hexdigest()
        assert actual == expected, (name, actual, expected)
    files = []
    for path in sorted(HERE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.name == 'manifest.json':
            continue
        if path.name.startswith('pilot'):
            continue
        data = path.read_bytes()
        files.append(dict(path=path.relative_to(HERE).as_posix(), bytes=len(data),
                          sha256=hashlib.sha256(data).hexdigest()))
    write('manifest.json', dict(environment=versions,
        scope='Independent reproducibility closure for the reviewer-revision forced stability and clearance scan.',
        inputs_verified_against_protocol=True, files=files,
        excluded=['__pycache__ (regenerable Numba/Python caches)',
                  'pilot* (preliminary performance/range check, not used in formal results)'],
        scientific_file_count=len(files), total_bytes=sum(r['bytes'] for r in files)))
    print(json.dumps(dict(files=len(files), total_bytes=sum(r['bytes'] for r in files))),flush=True)


if __name__ == '__main__': main()
