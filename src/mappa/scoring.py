"""The sole MAPPA TS scorer: significant asynchronous event-subset motion.

No overlap scorer, policy selector or fallback. The 0.4 review cue is not
implemented here: the score remains a descriptor, not a chemical acceptance gate.
"""
import numpy as np
from ase.data import atomic_numbers, covalent_radii
POLICY = 'significant_subset_motion_v1'
ACTIVE_SETTING = (1.5, 6)
ACTIVE_RECIPE = 'contact_1.5_6'
CONTACT_SETTINGS = [ACTIVE_SETTING]

def activity_balance(signed_derivatives):
    """Amplitude-derived event weights; inactive events carry zero weight.

    Squared positive bond-length derivatives give dimensionless normalized
    shares, not physical energies. Their inverse participation ratio is the
    effective number of co-active intended events (between 1 and m).
    """
    a = np.maximum(np.asarray(signed_derivatives, float), 0.)
    largest = float(np.max(a, initial=0.))
    if largest < 1e-12:
        return dict(event_activity_weights=np.zeros_like(a).tolist(),
                    effective_event_count=0., synchrony=0.)
    weights = (a/largest)**2
    weights /= weights.sum()
    effective = float(1./(weights@weights))
    synchrony = float(np.clip((effective-1)/(len(a)-1), 0, 1)) if len(a)>1 else 1.
    return dict(event_activity_weights=weights.tolist(), effective_event_count=effective,
                synchrony=synchrony)


def pair_geometry(xyz,elements):
    xyz=np.asarray(xyz,float);i,j=np.triu_indices(len(xyz),1)
    if not np.isfinite(xyz).all():raise ValueError('Nonfinite geometry')
    delta=xyz[j]-xyz[i];distance=np.linalg.norm(delta,axis=1)
    if np.any(distance<1e-8):raise ValueError('Coincident atoms')
    radii=np.array([covalent_radii[atomic_numbers[e]] for e in elements])
    lookup=np.full((len(xyz),len(xyz)),-1,int);lookup[i,j]=lookup[j,i]=np.arange(len(i))
    return dict(i=i,j=j,distance=distance,axis=delta/distance[:,None],
                radii_sum=radii[i]+radii[j],lookup=lookup)


def prepare_contacts(geometry,modes,settings=CONTACT_SETTINGS):
    out={}
    for k,mode in modes.items():
        q=np.array(mode['displacements'],float)
        if not np.isfinite(q).all() or np.linalg.norm(q)==0:raise ValueError('Finite nonzero mode required')
        q/=np.linalg.norm(q)
        velocity=np.einsum('ij,ij->i',geometry['axis'],q[geometry['j']]-q[geometry['i']])
        row={}
        for factor,power in settings:
            x=geometry['distance']/(factor*geometry['radii_sum']);xn=x**power
            gradient=power*xn/(geometry['distance']*(1+xn)**2)
            rate=gradient*velocity
            row[f'contact_{factor}_{power}']=dict(rate=rate,total=float(rate@rate),gradient=gradient)
        out[int(k)]=row
    return out


def contact_score(geometry,contacts,events,*,orientation=None):
    ix=np.array([geometry['lookup'][tuple(e['target_bond'])] for e in events]);sign=np.array([1 if e['kind']=='broken' else -1 for e in events])
    assert len(set(ix))==len(ix),'Duplicate event coordinates would double-count activity'
    values=contacts['rate'][ix]*sign;denominator=contacts['total'];choices=[]
    if orientation not in (None,1,-1):raise ValueError('Mode orientation must be +1 or -1')
    for sign_choice in ([1,-1] if orientation is None else [orientation]):
        a=sign_choice*values;positive=np.maximum(a,0);activity=float(positive@positive)
        balance=activity_balance(a);factor=.95+.05*balance['synchrony']
        base=float(np.sqrt(activity/denominator)) if denominator>1e-24 else 0.
        choices.append(dict(score=base*factor,base_score=base,reactive_amplitude=np.sqrt(activity),
            total_contact_amplitude=np.sqrt(denominator),orientation=sign_choice,
            signed_contact_rates=a.tolist(),synchrony_factor=factor,**balance))
    return max(choices,key=lambda x:(x['score'],x['orientation']))



def restrict_contacts(geometry, contacts, pairs):
    indices = np.array([geometry['lookup'][a, b] for a, b in pairs], int)
    return dict(contacts, total=float(np.sum(contacts['rate'][indices] ** 2)))

def evaluate_subset(geometry,contacts,events,displacement,*,reference_factor=1.):
    if reference_factor<=0:raise ValueError('Positive reference factor required')
    q=np.asarray(displacement,float)
    norm=np.linalg.norm(q)
    if not np.isfinite(q).all() or norm==0:raise ValueError('Finite nonzero displacement required')
    q=q/norm
    speed2=np.sum(q*q,axis=1)
    reference=float(np.sqrt(np.sum(speed2**2)/np.sum(speed2)))
    ix=np.array([geometry['lookup'][tuple(e['target_bond'])] for e in events])
    signs=np.array([1 if e['kind']=='broken' else -1 for e in events])
    derivative=signs*np.einsum('ij,ij->i',geometry['axis'][ix],q[geometry['j'][ix]]-q[geometry['i'][ix]])
    x=geometry['distance'][ix]/(1.5*geometry['radii_sum'][ix])
    proximity=np.minimum(1.,2./(1.+x**6))
    choices=[]
    for orientation in [1,-1]:
        d=contact_score(geometry,contacts,events,orientation=orientation)
        speeds=np.maximum(orientation*derivative,0.)*proximity
        support=speeds/(reference_factor*reference)
        strongest=int(np.argmax(support))
        strength=float(min(1.,support[strongest]))
        choices.append(dict(d,contact_alignment=d['base_score'],contact_ratio_score=d['score'],
            subset_strength=strength,reference_atomic_speed=reference,reference_factor=reference_factor,
            strongest_event=strongest,event_support_ratios=support.tolist(),event_proximity_factors=proximity.tolist(),
            event_supported_speeds=speeds.tolist(),base_score=strength*d['base_score'],score=strength*d['score']))
    result=max(choices,key=lambda d:(d['score'],d['orientation']))
    # Sensitivity at fixed assignment, with both common orientations reconsidered.
    result['reference_sensitivity_scores']={str(factor):max(
        min(1.,max(c['event_supported_speeds'])/(factor*reference))*c['contact_ratio_score']
        for c in choices) for factor in [.75,1.,1.25]}
    assert 0<=result['score']<=result['subset_strength']+1e-12<=1+2e-12
    return result


def subset_score(xyz,displacement,elements,events,pairs,*,reference_factor=1.):
    if not events:
        raise ValueError('At least one reaction event is required')
    geometry=pair_geometry(xyz,elements)
    contacts=prepare_contacts(geometry,{0:dict(displacements=displacement)},[ACTIVE_SETTING])[0][ACTIVE_RECIPE]
    return evaluate_subset(geometry,restrict_contacts(geometry,contacts,pairs),events,displacement,
        reference_factor=reference_factor)
