"""Production input, coefficient lifecycle, rollback and HDF metadata checks."""
from pathlib import Path
import subprocess
import sys
import tempfile
import h5py
import numpy as np

def main():
    executable=Path(sys.argv[1]).resolve()
    valid=[('omitted','',1.),('one','transport_1d_relaxation=1',1.),
           ('quarter','transport_1d_relaxation=.25',.25),('half','transport_1d_relaxation=.5',.5),
           ('small','transport_1d_relaxation=1e-6',1e-6)]
    invalid=['0','-.1','1.01','NaN','Infinity']
    message='transport_1d_relaxation must be finite and in (0,1]'
    with tempfile.TemporaryDirectory(prefix='transport-relaxation-') as root:
        for name,entry,expected in valid:
            d=Path(root)/name;d.mkdir()
            (d/'transport_model.nml').write_text('&TRANSPORT_MODEL_1D_LST\n'+entry+'\n/\n')
            r=subprocess.run([str(executable)],cwd=d,capture_output=True,text=True,timeout=30)
            if r.returncode or 'transport relaxation checks: PASS' not in r.stdout:
                raise AssertionError(f'{name}\n{r.stdout}{r.stderr}')
            with h5py.File(d/'transport_relaxation.h5') as h:
                assert np.asarray(h['transport_1d/params/transport_1d_relaxation']).item()==expected
                defect=h['transport_1d/diagnostics/coefficient_relative_defect'][()]
                assert defect.shape==(5,) and np.isfinite(defect).all() and (defect>=0).all()
                if expected==1:np.testing.assert_array_equal(defect,np.zeros(5))
                else:assert (defect>0).any()
                order=np.asarray(h['transport_1d/diagnostics/coefficient_order']).ravel()[0].decode().strip()
                assert order=='D,nu_mom,chi_i,chi_e,vpinch'
        for value in invalid:
            d=Path(root)/('invalid'+value);d.mkdir()
            (d/'transport_model.nml').write_text('&TRANSPORT_MODEL_1D_LST\ntransport_1d_relaxation='+value+'\n/\n')
            for args in [[],[value]]:
                r=subprocess.run([str(executable),*args],cwd=d,capture_output=True,text=True,timeout=30)
                assert r.returncode and message in r.stdout+r.stderr,(value,args,r.stdout,r.stderr)
    print(f'transport relaxation: {len(valid)+2*len(invalid)} input/lifecycle/metadata cases PASS')
if __name__=='__main__':main()
