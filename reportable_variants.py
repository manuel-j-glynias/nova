


def add_reason_to_variant(variant,reason):
    if not 'reasons' in variant:
        variant['reasons'] = []
    variant['reasons'].append(reason)

def add_lack_of_reason_to_variant(variant,reason):
    if not 'lack_of_reasons' in variant:
        variant['lack_of_reasons'] = []
    variant['lack_of_reasons'].append(reason)



def is_protein_altering(variant):
    if variant is None or 'is_protein_altering' not in variant:
        print('yikes')
    altering_ = variant['is_protein_altering']
    if not altering_:
        add_lack_of_reason_to_variant(variant, 'not a protein altering mutation')
    return altering_


def is_in_clinvar(variant):
    if variant['in_clinvar']:
        add_reason_to_variant(variant,'in ClinVar')
    else:
        add_reason_to_variant(variant, 'not in ClinVar')
    return variant['in_clinvar']


def is_clinvar_benign(variant):
    benign_ = variant['in_clinvar'] and variant['is_clinvar_benign']
    if benign_:
        add_reason_to_variant(variant, 'marked as benign in ClinVar')
    # else:
    #     add_reason_to_variant(variant, 'not marked as benign in ClinVar')
    return benign_


def is_clinvar_pathogenic(variant):
    pathogenic_ = variant['in_clinvar'] and variant['is_clinvar_pathogenic']
    if pathogenic_:
        add_reason_to_variant(variant, 'marked as pathogenic in ClinVar')
    else:
        add_lack_of_reason_to_variant(variant, 'not marked as pathogenic in ClinVar')
    return pathogenic_


def is_oncogene(variant):
    oncogene_ = variant['gene_category'] == 'Oncogene'
    if oncogene_:
        add_reason_to_variant(variant, 'in an oncogene')
    else:
        add_reason_to_variant(variant, 'in a tumor suppressor gene')
    return oncogene_

def is_ckb_gof(variant):
    gof = variant['is_gain_of_function']
    if gof:
        add_reason_to_variant(variant, 'annotated as gain of function in the JAX CKB')
    else:
        add_lack_of_reason_to_variant(variant, 'not annotated as gain of function in the JAX CKB')

    return gof


def at_hot_spot_position(variant):
    hotspots_ = 'hotspots' in variant and len(variant['hotspots']) > 0
    if hotspots_ :
        add_reason_to_variant(variant, 'annotated as a hot spot variant')
    else:
        add_lack_of_reason_to_variant(variant, 'not annotated as a hot spot variant')

    return hotspots_


def is_ckb_lof(variant):
    lof = variant['is_loss_of_function']
    if lof:
        if variant['is_loss_of_function']:
            add_reason_to_variant(variant, 'annotated as loss of function in the JAX CKB')
        else:
            add_lack_of_reason_to_variant(variant, 'not annotated as loss of function in the JAX CKB')
    return lof


def is_truncating(variant):
    truncating = variant['is_truncating_variants']
    if truncating:
        add_reason_to_variant(variant, 'is a truncating mutation')
    else:
        add_lack_of_reason_to_variant(variant, 'is not a truncating mutation')

    return truncating


def is_sift_deleterious_and_polyphen_damaging(variant):
    deleterious_ = variant['predicted_deleterious']
    if deleterious_:
        add_reason_to_variant(variant, 'predicted deleterious by SIFT (' + variant['sift_prediction'] + ') and ' \
                                'PolyPhen-2 (' + variant['polyphen_prediction'] + ')')
    else:
        add_lack_of_reason_to_variant(variant, 'not predicted deleterious by SIFT and PolyPhen-2')

    return deleterious_


def is_near_gof_or_lof_mutations(variant):
    near = variant['is_near_LOF_mutation']
    if near:
        add_reason_to_variant(variant, 'in a region containing known deleterious mutations')
    else:
        add_lack_of_reason_to_variant(variant, 'not in a region containing known deleterious mutations')
    return near

def mark_as_inactivating(variant):
    variant['is_oncogene'] = False
    variant['is_activating'] = False
    variant['is_inactivating'] = True
    variant['report_status'] = 'reportable'


def mark_as_activating(variant):
    variant['is_oncogene'] = True
    variant['is_activating'] = True
    variant['is_inactivating'] = False
    variant['report_status'] = 'reportable'

def is_cnv_gain(variant):
    return variant['pdot'] == 'amp'

def is_cnv_loss(variant):
    return variant['pdot'] == 'loss'

def is_cnv_reportable(variant):
    if is_oncogene(variant):
        if is_cnv_gain(variant):
            mark_as_activating(variant)
    else:
        if is_cnv_loss(variant):
            mark_as_inactivating(variant)


def not_clinvar_pathogenic(variant):
    if is_oncogene(variant):
        if is_ckb_gof(variant) or is_ckb_lof(variant):
            mark_as_activating(variant)
        else:
            if at_hot_spot_position(variant):
                mark_as_activating(variant)
            else:
                mark_as_not_reported(variant)
    else:  # is tumor suppressor
        if is_ckb_gof(variant) or is_ckb_lof(variant):
            mark_as_inactivating(variant)
        else:
            if is_truncating(variant):
                mark_as_inactivating(variant)
            else:
                if is_sift_deleterious_and_polyphen_damaging(variant):
                    mark_as_inactivating(variant)
                else:
                    if is_near_gof_or_lof_mutations(variant):
                        mark_as_inactivating(variant)
                    else:
                        mark_as_not_reported(variant)


def is_snv_reportable(variant):
    if is_protein_altering(variant):
        if is_in_clinvar(variant):
            if is_clinvar_benign(variant):
                mark_as_not_reported(variant)
            else:
                if is_clinvar_pathogenic(variant):
                    variant['report_status'] = 'reportable'
                    if is_oncogene(variant):
                        mark_as_activating(variant)
                    else:
                        mark_as_inactivating(variant)
                else:   #not clinvar pathogenic
                    not_clinvar_pathogenic(variant)
        else:
            not_clinvar_pathogenic(variant)
    else:
        mark_as_not_reported(variant)

def mark_as_not_reported(variant):
    variant['report_status'] = 'not_reported'
    variant['category'] = 'VUS'



