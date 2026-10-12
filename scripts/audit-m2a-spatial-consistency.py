"""Continuum Stokes fields through current J and stress sampling; no CFD solve."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import surface_reference as ref

spec=importlib.util.spec_from_file_location('wall_reference',Path(__file__).with_name('audit-m2a-wall-operators.py'))
wall=importlib.util.module_from_spec(spec);spec.loader.exec_module(wall)


def fixture(mode,n,phase):
    h=.12/n;radius=.01;centre=[phase*h]*3;origin=[-.06+h/2]*3
    count=math.ceil(4*math.pi*radius**2/(2.25*h*h));errors=[]
    omega=[0.,0.,.1 if mode=='rotate' else 0.]
    def field(index):
        r=[o+i*h-c for o,i,c in zip(origin,index,centre)]
        return ref.cross(omega,r)+[0.] if sum(x*x for x in r)<=radius**2 else wall.analytic(r,mode)
    for i in range(count):
        z=1-2*(i+.5)/count;phi=math.pi*(3-math.sqrt(5))*i
        normal=[math.sqrt(1-z*z)*math.cos(phi),math.sqrt(1-z*z)*math.sin(phi),z]
        point=[c+radius*x for c,x in zip(centre,normal)]
        actual=wall.sample(point,origin,h,field)[:3]
        target=ref.cross(omega,[radius*x for x in normal])
        errors.append(math.sqrt(sum((a-b)**2 for a,b in zip(actual,target))))
    stress=wall.fixture(mode,n,phase,'linear')
    return {'markers':count,'wall_rms_m_s':math.sqrt(sum(x*x for x in errors)/count),
            'wall_max_m_s':max(errors),'stress_load_ratio':stress['stress_load_ratio'],
            'reference_speed_m_s':.01 if mode=='fixed' else .001}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data={'scope':__doc__,'cases':[{'mode':mode,'mesh_n':n,'phase_h':phase,**fixture(mode,n,phase)}
         for mode in ('fixed','rotate') for n in (36,48,64) for phase in (0.,.25)]}
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
