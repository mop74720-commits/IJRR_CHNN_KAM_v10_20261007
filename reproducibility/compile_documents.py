"""Compile the supplement before the main article to resolve external references."""
from pathlib import Path
import subprocess,re

root=Path(__file__).resolve().parent.parent
for name in ('supplement','main'):
    for run in range(3):
        result=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',
            '-file-line-error','-no-shell-escape',name+'.tex'],cwd=root,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if result.returncode:
            print(result.stdout.decode('utf-8',errors='replace')[-6000:])
            raise SystemExit(result.returncode)
        log=(root/(name+'.log')).read_text(encoding='utf-8',errors='replace')
        if run>=1 and not re.search(r'Rerun to get|Label\(s\) may have changed|There were undefined',log):
            break
    warnings=re.findall(r'^.*(?:Overfull|undefined|LaTeX Error|Rerun to get|Label\(s\) may have changed).*',log,re.M)
    if warnings:
        print('\n'.join(warnings));raise SystemExit(1)
    print(name+': '+re.search(r'Output written on .+?\.pdf \([^\n]+',log)[0])
