

def get_pertinent_negatives_for_report(patient,variant_groups_dict):
    reportable_disease = patient['reportable_disease_name']
    if reportable_disease in variant_groups_dict:
        entries = variant_groups_dict[reportable_disease]
    else:
        entries = []
    entries.extend(variant_groups_dict['All Solid Tumor'])
    negatives = []
    v2g_dict = {}
    for entry in entries:
        vg = entry['variant_group']
        v2g_dict[entry['variant_name']]= entry['variant_group']
        if not vg in negatives:
            negatives.append(vg)

    for variant in patient['ckb_variants']:
        if variant in v2g_dict:
            vg = v2g_dict[variant]
            if vg in negatives:
                negatives.remove(vg)

    negatives.sort()
    return negatives