"""Check the Lean pilot, its input provenance, and rejection of a false bound."""

# The source mutation matches Lean's rational type notation literally.
# ruff: noqa: RUF001

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ALLOWED_AXIOMS = {'propext', 'Classical.choice', 'Quot.sound'}


def run(command):
    return subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lake', default='lake', help='Lake executable (default: on PATH)')
    args = parser.parse_args()
    print(run([sys.executable, 'export_certificate.py', '--check']).stdout, end='')
    build = run([args.lake, 'build'])
    print(build.stdout, end='')
    print(build.stderr, end='', file=sys.stderr)
    audit = run([args.lake, 'env', 'lean', 'Audit.lean'])
    print(audit.stdout, end='')
    axioms = re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", audit.stdout)
    expected = re.findall(r'^#print axioms (\S+)', (ROOT / 'Audit.lean').read_text(), re.MULTILINE)
    if {name for name, _ in axioms} != set(expected):
        raise RuntimeError('The axiom audit did not report every requested theorem')
    for name, dependencies in axioms:
        used = {item.strip() for item in dependencies.split(',') if item.strip()}
        if used - ALLOWED_AXIOMS:
            raise RuntimeError(f'{name} uses unexpected axioms: {used - ALLOWED_AXIOMS}')

    source = (ROOT / 'StepDetection' / 'SavedCertificate.lean').read_text()
    changed, count = re.subn(
        r'^def displayUpper : ℚ := .*$', 'def displayUpper : ℚ := 0', source, flags=re.MULTILINE
    )
    if count != 1:
        raise RuntimeError('Could not locate the bound to corrupt')
    # Keep the false proof outside the library and never overwrite the saved source.
    with tempfile.TemporaryDirectory(prefix='negative-', dir=ROOT / '.lake') as directory:
        path = Path(directory) / 'RejectedCertificate.lean'
        path.write_text(changed)
        rejected = subprocess.run(
            [args.lake, 'env', 'lean', str(path)], cwd=ROOT, text=True, capture_output=True
        )
        if rejected.returncode == 0 or 'unsolved goals' not in rejected.stdout:
            raise RuntimeError(
                f'Expected an unprovable arithmetic claim: {rejected.stdout}{rejected.stderr}'
            )
    print('Corrupted certificate (upper bound changed to zero): rejected by Lean.')
    print('Lean pilot verification passed.')


if __name__ == '__main__':
    try:
        main()
    except subprocess.CalledProcessError as error:
        print(error.stdout, end='', file=sys.stderr)
        print(error.stderr, end='', file=sys.stderr)
        raise SystemExit(error.returncode) from error
