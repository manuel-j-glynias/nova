import utils
import read_variants
import annotate
import reportable_variants
import pprint
import report_generator
import evidence
import os
import shutil
import negatives
import pathologist_comments
import additional_io



def get_variant_file_path():
    path = ''
    files = utils.get_list_of_files('input')
    for f in files:
        if 'VariantSummaries' in f and f.endswith('.csv'):
            path = f
            break
    # path = 'input/VariantSummaries_OA_2019_11_04.csv'
    return path


def get_immune_results_file_path():
    # path = 'input/ImmuneResultsSummaries_OA_2019_11_04.csv'
    path = ''
    files = utils.get_list_of_files('input')
    for f in files:
        if 'ImmuneResultsSummaries' in f and f.endswith('.csv'):
            path = f
            break
    return path


def find_all_reportable_markers(patient):
    patient['reportable_markers'] = []
    for v in patient['snv']:
        v['type'] = 'snv'
        if v['report_status'] == 'reportable':
            patient['reportable_markers'].append(v)
    for v in patient['cnv']:
        v['type'] = 'cnv'
        if v['report_status'] == 'reportable':
            patient['reportable_markers'].append(v)
    for v in patient['io']:
        if v['report_status'] == 'reportable':
            patient['reportable_markers'].append(v)


def variant_name(variant):
    return variant['gene'] + ' ' + variant['pdot']


def get_markers(patient):
    indications = patient['indications']
    indications.extend(patient['off_label'])
    indications.extend(patient['contraindications'])
    for indication in indications:
        if indication['approval_index'] == 1:
            for variant in indication['variants']:
                if variant['type'] != 'io' and variant['pdot'] != 'wild-type' \
                        and variant_name(variant) not in patient['cdx_genomic_markers']:
                    patient['cdx_genomic_markers'].append(variant_name(variant))
                elif variant['type'] == 'io' and variant_name(variant) not in patient['cdx_io']:
                        patient['cdx_io'].append(variant_name(variant))
        elif indication['approval_index'] <= 5:
            for variant in indication['variants']:
                if variant['type'] != 'io' and variant['pdot'] != 'wild-type' \
                        and variant_name(variant) not in patient['cdx_genomic_markers'] \
                        and variant_name(variant) not in patient['actionable_genomic_markers']:
                    patient['actionable_genomic_markers'].append(variant_name(variant))


def variants_name(variant_list):
    name = ''
    for variant in variant_list:
        v_name = variant_name(variant)
        if len(name) > 0:
            name = name + '+'
        name = name + v_name
    return name

def append_unique_therapy(therapy_list,therapy):
    found = False
    for t in therapy_list:
        if t['therapy']==therapy['therapy']:
            found = True
            break
    if not found:
        therapy_list.append(therapy)

def get_therapies_by_marker(indications):
    therapies_by_marker = {}
    convert_approval_index = ['','CDx', 'Contraindicated', 'NCCN','Phase III','Phase II']
    for indication in indications:
        v_name = variants_name(indication['variants'])

        if v_name not in therapies_by_marker:
            therapies_by_marker[v_name] = {'variants_name':v_name, 'therapies':[]}
        variant_entry = therapies_by_marker[v_name]
        therapy = {'therapy': indication['therapy']}
        therapy['setting'] = '-'
        if indication['contraindicated']:
            therapy['setting'] = 'Contraindicated'
        elif indication['off_label']:
            therapy['setting'] = 'Off Label'
        therapy['evidence'] = convert_approval_index[indication['approval_index']] + ' (' + indication[
            'ampCapAscoInferredTier'] + indication['ampCapAscoEvidenceLevel'] + ')'
        append_unique_therapy(variant_entry['therapies'],therapy)
    list_of_markers_with_therapies = []
    for key, value in therapies_by_marker.items():
        list_of_markers_with_therapies.append({'marker':key, 'therapies':value['therapies']})
    return list_of_markers_with_therapies


def front_page_data(patient,evidence_collection):
    patient['cdx_genomic_markers'] = []
    patient['actionable_genomic_markers'] = []
    patient['pertinent_negatives'] = negatives.get_pertinent_negatives_for_report(patient,evidence_collection)
    get_markers(patient)
    indications = patient['contraindications']
    indications.extend(patient['indications'])
    indications.extend(patient['off_label'])
    patient['therapies_by_marker'] = get_therapies_by_marker(indications)
    # if len(patient['therapies_by_marker'])==0:
    #     patient['therapies_by_marker'] = get_therapies_by_marker(patient['off_label'])


def third_page_data(patient):
    therapy_details_dict = {}
    indications = patient['indications']
    indications.extend(patient['off_label'])
    indications.extend(patient['contraindications'])
    for indication in indications:
        for variant in indication['variants']:
            if variant['pdot'] != 'wild-type':
                name = variant_name(variant)
                if name not in therapy_details_dict:
                    therapy_details_dict[name] = {'variant':variant, 'strong':[], 'moderate':[], 'emerging':[] }
                therapy_detail = therapy_details_dict[name]
                if indication['off_label']:
                    indication['availability'] = 'Off Label'
                cap = indication['ampCapAscoInferredTier'] + indication['ampCapAscoEvidenceLevel']
                if cap == 'IA':
                    append_unique_therapy(therapy_detail['strong'],indication)
                elif cap == 'IB':
                    append_unique_therapy(therapy_detail['moderate'],indication)
                else:
                    append_unique_therapy(therapy_detail['emerging'],indication)

    therapy_details = therapy_details_dict.values()
    # for t_by_m in patient['therapies_by_marker']:
    #     print(t_by_m)
    #     therapy_detail = {'variant':'EGFR L858R (c.2573T>G)', 'gene_description':'Description', 'vaf':'19.7'}
    #     therapy_details.append(therapy_detail)
    patient['therapy_details'] = therapy_details

def handle_one_patient(patient, db, strands):
    patient['snv'] = []
    patient['cnv'] = []
    patient['io'] = []
    read_variants.add_io_data(patient)
    for variant in patient['unannotated_snv']:
        annotated_variant = annotate.annotate_snv(strands, variant, db)
        reportable_variants.is_snv_reportable(annotated_variant)
        patient['snv'].append(annotated_variant)
    for variant in patient['unannotated_cnv']:
        annotated_variant = annotate.annotate_cnv(variant, db)
        reportable_variants.is_cnv_reportable(annotated_variant)
        patient['cnv'].append(annotated_variant)
    find_all_reportable_markers(patient)
    evidence.find_therapies(patient, db['evidence'])
    front_page_data(patient, db['evidence'])
    third_page_data(patient)
    additional_io.add_additional_io(patient)
    # pprint.pprint(patient)

def create_one_report(patient):
    report_generator.render_html_report(patient)



def create_recommendations(patient, db):
    patient['pathologist_comments'] = pathologist_comments.get_draft_comments(patient)


def empty_output_dir():
    for root, dirs, files in os.walk('output'):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))

def main():
    empty_output_dir()
    client = utils.get_mongo_client()
    db = utils.get_database(client,'omni')
    patients = read_variants.read_immune_results_file(get_immune_results_file_path())
    read_variants.read_all_variants(patients, get_variant_file_path())
    strands = annotate.read_strands('data/strands.xlsx')

    for order_id in patients.keys():
        patient = patients[order_id]
        read_variants.add_patient_data(patient)
        print(patient['fake_order_id'], '=',order_id)
        handle_one_patient(patient, db, strands)
        create_recommendations(patient,db)
        create_one_report(patient)
        # break
        




if __name__ == "__main__":
    main()
