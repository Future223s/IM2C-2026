"""Predict category percentages from a JSON environmental vector."""
import csv,json,sys
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'results/tree_ridge_coefficients.csv'
def predict(values):
    with p.open(encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f); cats=reader.fieldnames[3:]; rows=list(reader)
    unknown=set(values)-{r['Predictor'] for r in rows[1:]}
    if unknown: raise ValueError(f'Unknown predictors: {unknown}')
    result={c:float(rows[0][c]) for c in cats}
    for r in rows[1:]:
        x=values.get(r['Predictor']); z=0 if x is None else (float(x)-float(r['Mean']))/float(r['SD'])
        for c in cats: result[c]+=z*float(r[c])
    total=sum(max(0,v) for v in result.values())
    return {'affine_percent':result,'normalized_percent':{c:100*max(0,v)/total for c,v in result.items()}}
if __name__=='__main__':
    print(json.dumps(predict(json.loads(Path(sys.argv[1]).read_text())),indent=2))
