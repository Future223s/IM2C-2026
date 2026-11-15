"""Rebuild FIA distributions and lambda=3 ridge outputs. Run from any directory."""
from pathlib import Path
import csv
import numpy as np
import openpyxl
from model_distribution import distribution
ROOT=Path(__file__).resolve().parents[1]
CITIES=['Montpelier VT','Old Forge NY','Binghamton NY','Scranton PA','State College PA','Frostburg MD','Valley Head WV','Shenandoah VA','Roanoke VA','Boone NC','Asheville NC']
FILES=['montpelier_vt','oldforge_ny','binghamton_ny','scranton_pa','statecollege_pa','frostburg_md','valleyhead_wv','shenandoah_va','roanoke_va','boone_nc','asheville_nc']
def write(path,header,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)
def fit(x,y):
    mu=np.nanmean(x,axis=0); x=np.where(np.isnan(x),mu,x); sd=x.std(axis=0); sd=np.where(sd==0,1,sd); z=(x-mu)/sd
    b=y.mean(axis=0); B=np.linalg.solve(z.T@z+3*np.eye(z.shape[1]),z.T@(y-b));return mu,sd,b,B
def predict(x,model):
    mu,sd,b,B=model;return (np.where(np.isnan(x),mu,x)-mu)/sd@B+b
def main():
    out=ROOT/'results'; data=ROOT/'data/processed'; rows=[]
    for city,name in zip(CITIES,FILES):
        path=ROOT/'data/raw_fia'/(name+'_data.csv' if name!='asheville_nc' else 'asheville_nc_data(1).csv')
        d=distribution(path); rows.extend((city,*v) for v in d)
    write(data/'category_counts.csv',['City','Category','Estimated trees','Percentage'],rows)
    s=list(openpyxl.load_workbook(data/'data_environment.xlsx',data_only=True).worksheets[0].values)
    cats=[r[0] for r in s[1:12]]; y=np.array([r[1:] for r in s[1:12]],float).T
    features=[r for r in s[13:] if r[0] not in ['Town elevation (m)','Autumn equinox daylight (h)']]; names=[r[0] for r in features]; x=np.array([r[1:] for r in features],float).T
    actual=np.array([r[3] for r in rows]).reshape(11,11); np.testing.assert_allclose(actual,y,atol=1e-7)
    model=fit(x,y); mu,sd,b,B=model; fitted=predict(x,model); cv=np.empty_like(y)
    for i in range(len(y)):
        m=np.arange(len(y))!=i; cv[i]=predict(x[i:i+1],fit(x[m],y[m]))[0]
    write(out/'recomputed_coefficients.csv',['Predictor','Mean','SD',*cats],[['Intercept','','',*b]]+[[n,m,d,*v] for n,m,d,v in zip(names,mu,sd,B)])
    write(out/'recomputed_predictions.csv',['City','Category','Actual','Fitted','LOOCV'],[(city,cats[j],y[i,j],fitted[i,j],cv[i,j]) for i,city in enumerate(CITIES) for j in range(len(cats))])
    sst=((y-y.mean(0))**2).sum(0); tr=1-((y-fitted)**2).sum(0)/sst; cr=1-((y-cv)**2).sum(0)/sst; rm=np.sqrt(((y-cv)**2).mean(0))
    write(out/'ridge_metrics.csv',['Category','Training R2','LOOCV R2','LOOCV RMSE (percentage points)'],zip(cats,tr,cr,rm))
    print('Verified raw FIA distributions against workbook; generated ridge results.')
if __name__=='__main__':main()
