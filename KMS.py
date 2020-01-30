import requests
import urllib.parse
import pprint
import urllib3


def call_go_api(ext):
    val = None
    server = "https://10.10.200.107/api/"
    server_ext = server + ext
    r = requests.get(server_ext, headers={"Content-Type": "application/json", 'Authorization': 'Token 699822a95346d78da5c0e2c317279ae2a1269f32'}, verify=False)
    if not r.ok:
        #         r.raise_for_status()
        print(r)
    else:
        val = r.json()
    return val

def get_go_alteration_from_ckb_id(ckb_id):
    alteration_name = None
    q = urllib.parse.quote('JAX' + ckb_id)
    ext = 'alterations/?codes=JAX-CKB%3A%20' + q
    val = call_go_api(ext)
    if 'results' in val:
        if len(val['results'])>0:
            alteration_name = val['results'][0]['name']
    return alteration_name

def get_go_disease_from_ckb_id(ckb_id):
    disease_name = None
    q = urllib.parse.quote('JAX'+ckb_id)
    ext = 'diseases/?codes=JAX-CKB%3D' + q
    val = call_go_api(ext)
    if 'results' in val:
        if len(val['results'])>0:
            disease_name = val['results'][0]['name']
    return disease_name

def get_kms_disease_from_omni_disease(disease):
    disease_name = None
    q = urllib.parse.quote(disease)
    ext = 'diseases/suggest?q=' + q
    val = call_go_api(ext)
    if 'results' in val:
        if len(val['results'])>0:
            disease_name = val['results'][0]['text']
    return disease_name



def get_therapies(ext,alt_dict):
    therapies = []
    negatives = []
    val = call_go_api(ext)
    if 'results' in val:
        if len(val['results']) > 0:
            for therapy in val['results']:
                therapy_dict = {}
                therapy_dict['drugs'] = therapy['drugs']
                therapy_dict['responses'] = therapy['response_list']
                therapy_dict['response_is_positive'] = therapy['response_is_positive']
                if 'other_condition' in therapy:
                    therapy_dict['other'] = therapy['other_condition']
                if 'amp_tier_evidence' in therapy:
                    therapy_dict['amp_tier_evidence'] = therapy['amp_tier_evidence']
                therapy_dict['detected_alterations'] = []
                if 'setting_source_pairs' in therapy:
                    therapy_dict['setting_source_pairs'] = therapy['setting_source_pairs']
                if 'drug_categories' in therapy:
                    therapy_dict['drug_categories'] = therapy['drug_categories']
                if 'detected_alterations' in therapy:
                    for detected in therapy['detected_alterations']:
                        for trigger in detected['trigger_alterations']:
                            if trigger in alt_dict:
                                variant = alt_dict[trigger]
                                therapy_dict['detected_alterations'].append(variant)
                if 'evidence_references' in therapy:
                    evidences = []
                    for evidence in therapy['evidence_references']:
                        evidence_dict = {'evidence_category': '-',
                                         'reference_title': '-',
                                         'reference_source_name': '-',
                                         'reference_source_id':'-'}
                        if 'evidence_category' in evidence:
                            evidence_dict['evidence_category'] = evidence['evidence_category']
                        if 'reference_title' in evidence:
                            evidence_dict['reference_title'] = evidence['reference_title']
                        if 'reference_source_name' in evidence:
                            evidence_dict['reference_source_name'] = evidence['reference_source_name']

                        if 'reference_source_id' in evidence:
                            evidence_dict['reference_source_id'] = evidence['reference_source_id']

                        evidences.append(evidence_dict)
                    therapy_dict['evidences'] = evidences
                if therapy_dict['response_is_positive']:
                    therapies.append(therapy_dict)
                else:
                    negatives.append(therapy_dict)
                therapy_dict['FDA_Approved'] = False
                if 'sources' in therapy and 'FDA' in therapy['sources']:
                    therapy_dict['FDA_Approved'] = True
                    if 'amp_tier_evidence' not in therapy_dict:
                        therapy_dict['amp_tier_evidence'] = '1A'
                if 'amp_tier_evidence' not in therapy_dict:
                    therapy_dict['amp_tier_evidence'] = '1B'
    return therapies,negatives

def match_go_MCG_therapies(diseases,alterations,alt_dict):
    ext = 'therapies/matches?data_set=MCG&sources=FDA&sources=NCCN'
    for disease in diseases:
        disease_name = urllib.parse.quote_plus(disease)
        ext += '&diseases=' + disease_name
    for alteration in alterations:
        alteration_name = urllib.parse.quote_plus(alteration)
        ext += '&alterations=' + alteration_name
    therapies,negatives =get_therapies(ext,alt_dict)
    return therapies, negatives


def match_go_CKB_therapies(diseases,alterations,tier_list,alt_dict):
    ext = 'therapies/matches?data_set=JAX-CKB'
    if tier_list is not None:
        for tier in tier_list:
            ext += '&amp_tier_evidence=' + tier
    for disease in diseases:
        disease_name = urllib.parse.quote_plus(disease)
        ext += '&diseases=' + disease_name
    for alteration in alterations:
        alteration_name = urllib.parse.quote_plus(alteration)
        ext += '&alterations=' + alteration_name
    therapies,negatives = get_therapies(ext,alt_dict)
    return therapies,negatives



def match_go_local_trials(diseases,alterations,zip,clincal_trial_distance,dob,alt_dict):
    dob = dob.replace('/','%2F')
    ext = 'trials/matches?zip=' + zip + '&distance_miles=' + str(clincal_trial_distance) + '&date_of_birth=' + dob
    for disease in diseases:
        disease_name = urllib.parse.quote_plus(disease)
        ext += '&diseases=' + disease_name
    for alteration in alterations:
        alteration_name = urllib.parse.quote_plus(alteration)
        ext += '&alterations=' + alteration_name
    val = call_go_api(ext)
    trials = []
    if 'results' in val:
        if len(val['results'])>0:
            for trial in val['results']:
                trial_dict = {'nct_id': trial['nct_id'],
                            'title': trial['title'],'phase':trial['phase'],'drugs':[],
                              'trigger_alterations':[], 'locations':[], 'evidences':[]}
                match_results = trial['match_results'][0]
                if 'treatment_contexts' in match_results and len(match_results['treatment_contexts'])>0:
                    for drug_dict in match_results['treatment_contexts'][0]['drugs']:
                        trial_dict['drugs'].append(drug_dict['name'])
                trials.append(trial_dict)
                for alt in trial['detected_alterations']:
                    for a in alt['trigger_alterations']:
                        if a in alt_dict:
                            variant = alt_dict[a]
                            trial_dict['trigger_alterations'].append(variant)
                facilities = []
                for location in trial['trial_locations']:
                    if location['facility_name'] not in facilities:
                        loc = {'facility_name':location['facility_name'],'city':location['city'],'state':location['state']}
                        trial_dict['locations'].append(loc)
                        facilities.append(location['facility_name'])
    return trials

def match_go_local_trials_by_drug(drug,zip,clincal_trial_distance,dob):
#     trials/matches?diseases=ANY&skip_alteration_match=true&drugs=crizotinib&zip=14203&date_of_birth=09%2F01%2F1955&distance_miles=150
    dob = dob.replace('/','%2F')
    ext = 'trials/matches?diseases=ANY&skip_alteration_match=true&zip=' + zip + '&distance_miles=' + str(
        clincal_trial_distance) + '&date_of_birth=' + dob + '&drugs=' + drug
    val = call_go_api(ext)
    trials = []
    if 'results' in val:
        if len(val['results']) > 0:
            for trial in val['results']:
                trial_dict = {'nct_id': trial['nct_id'],
                              'title': trial['title'], 'phase': trial['phase'], 'drugs': [],
                              'trigger_alterations': [], 'locations': [], 'evidences': []}
                match_results = trial['match_results'][0]
                if 'treatment_contexts' in match_results and len(match_results['treatment_contexts']) > 0:
                    for drug_dict in match_results['treatment_contexts']:
                        for drug in drug_dict['drugs']:
                            if drug['name'] not in trial_dict['drugs']:
                                trial_dict['drugs'].append(drug['name'])
                trials.append(trial_dict)
                facilities = []
                for location in trial['trial_locations']:
                    if 'facility_name' in location:
                        if location['facility_name'] not in facilities:
                            loc = {'facility_name': location['facility_name'], 'city': location['city'],
                                   'state': location['state']}
                            trial_dict['locations'].append(loc)
                            facilities.append(location['facility_name'])
    return trials
def variant_name(variant):
    if 'reported_variant' in variant:
        name = variant['reported_variant']
    else:
        name = variant['gene'] + ' ' + variant['pdot']
    return name


def handle_io(v,alt_list,alt_dict):
    if v['gene']=='PD-L1':
        if v['pdot']=='positive':
            alt_list.append('PD-L1 Expression')
            alt_dict['PD-L1 Expression'] = v
            alt_list.append('PD-L1 High Expression')
            alt_dict['PD-L1 High Expression'] = v
    elif v['gene']=='TMB':
        if v['pdot']=='high':
            alt_list.append('TMB-High')
            alt_dict['TMB-High'] = v
    elif v['gene']=='MSI':
        if v['pdot']=='high':
            alt_list.append('MSI-High')
            alt_dict['MSI-High'] = v
        elif v['pdot']=='stable':
            alt_list.append('MSI-Low')
            alt_dict['MSI-Low'] = v



def get_kms_nonsense_variant(gene):
    go_alt = None
    q = urllib.parse.quote(gene)
    ext = 'alterations/suggest?q=' + q + '%20Nonsense'
    val = call_go_api(ext)
    if 'results' in val:
        if len(val['results'])>0:
            go_alt = gene + ' Nonsense'
    if go_alt==None:
        ext = 'alterations/suggest?q=' + q + '%20Mutation'
        val = call_go_api(ext)
        if 'results' in val:
            if len(val['results']) > 0:
                go_alt = gene + ' Mutation'
    return go_alt


def key_for_therapy(therapy):
    s = ''
    for d in therapy['drugs']:
        if len(s)>0:
            s += '+'
        s += d.lower()
    return s


def find_therapies_and_clinical_trials(patient):
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    disease_list = []
    if 'ckb_disease_ids' in patient and patient['ckb_disease_ids'] is not None:
        for disease_id in patient['ckb_disease_ids']:
            disease_name = get_go_disease_from_ckb_id(disease_id)
            if disease_name not in disease_list and disease_name is not None:
                disease_list.append(disease_name)
    kms_disease_name = get_kms_disease_from_omni_disease(patient['OmniDisease'])
    if kms_disease_name is not None and kms_disease_name not in disease_list:
        disease_list.append(kms_disease_name)

    alt_list = []
    alt_dict = {}
    for v in patient['reportable_markers']:
        name = variant_name(v)
        if name.startswith('PD-L1') or name.startswith('TMB') or name.startswith('MSI'):
            handle_io(v,alt_list,alt_dict)
        elif name.startswith('BRCA') or name.startswith('TP53'):
            # alt_list.append(v['gene'] + " Mutation")
            alt_list.append(v['gene'] + " Inactivating Mutation")
            alt_dict[v['gene'] + " Inactivating Mutation"] = v
        else:
            go_alt = None
            if v['pdot'].endswith('*') or  v['pdot'].endswith('X'):
                go_alt = get_kms_nonsense_variant(v['gene'])
            if go_alt==None:
                if 'ckb_id' in v and len(v['ckb_id'])>0:
                    go_alt = get_go_alteration_from_ckb_id(v['ckb_id'])
                else:
                    print('no ckb for:',name)
            if go_alt is not None and go_alt not in alt_list:
                alt_list.append(go_alt)
                alt_dict[go_alt] = v

    # print('disease_list=',disease_list)
    # print('alt_list=',alt_list)
    therapies,negatives = match_go_MCG_therapies(disease_list,alt_list,alt_dict)
    therapy_dict = {}
    # print('num MCG therapies =',len(therapies))
    for therapy in therapies:
        key = key_for_therapy(therapy)
        if key in therapy_dict:
            print('duplicate key for:',key)
        else:
            therapy_dict[key] = therapy
        # print(key)
        # pprint.pprint(therapy)
    # if len(negatives)>0:
    #     print('num MCG negatives =', len(negatives))
    #     pprint.pprint(negatives)
    patient['kms_therapies'] = therapies
    patient['kms_negatives'] = negatives

    zip = '14203'
    clincal_trial_distance = 150

    trials = match_go_local_trials(disease_list,alt_list,zip,clincal_trial_distance,patient['DOB'],alt_dict)
    patient['kms_trials'] = trials
    trials_dict = {}
    for trial in trials:
        for drug in trial['drugs']:
            key = drug.lower()
            if key in trials_dict:
                trials = trials_dict[key]
                trials.append(trial)
            else:
                trials_dict[key] = [trial]
    # print(trials_dict.keys())
    # print('num trials =',len(trials))
    # pprint.pprint(trials)

    ckb_therapies,ckb_negatives = match_go_CKB_therapies(disease_list,alt_list,['1A','1B','2C','2D'],alt_dict)
    # print('num ckb therapies =',len(ckb_therapies))
    for ckb_therapy in ckb_therapies:
        if len(ckb_therapy['evidences'])>0:
            key = key_for_therapy(ckb_therapy)
            if key in therapy_dict:
                therapy = therapy_dict[key]
                if therapy['FDA_Approved']:
                    get_evdience_for_fda(ckb_therapy, therapy)
                else:
                    for evidence in ckb_therapy['evidences']:
                        therapy['evidences'].append(evidence)
                if 'other' in ckb_therapy:
                    if not 'other' in therapy:
                    #     therapy['other'] = therapy['other'] +  ' ' + ckb_therapy['other']
                    # else:
                        therapy['other'] = ckb_therapy['other']
            else:
                # print(key)
                if key in trials_dict:
                    # print(ckb_therapy)
                    for trial in trials_dict[key]:
                        for evidence in ckb_therapy['evidences']:
                            trial['evidences'].append(evidence)
    # pprint.pprint(therapies)
    # if len(therapies)==0 and len(ckb_therapies)>0:
    #     pprint.pprint(ckb_therapies)
    # if len(ckb_negatives)>0:
    #     print('num CKB negatives =', len(ckb_negatives))
    #     pprint.pprint(ckb_negatives)


def get_evdience_for_fda(ckb_therapy, therapy):
    for evidence in ckb_therapy['evidences']:
        if evidence['evidence_category'].startswith('FDA') and evidence['reference_source_name'] == 'PubMed':
            therapy['evidences'].append(evidence)
    if len(therapy['evidences'])==0:
        for evidence in ckb_therapy['evidences']:
            therapy['evidences'].append(evidence)



