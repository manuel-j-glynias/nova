import numbers
import csv
from random import randint
from faker import Faker

def read_fusion(row):
    variant = {}
    variant['HGNC_Symbol'] = row['HGNC_Symbol']
    variant['gene'] = row['HGNC_Symbol']
    variant['Variant'] = row['Variant']
    variant['reported_variant'] = variant['HGNC_Symbol'] + ' ' +variant['Variant']
    variant['pdot'] = 'rearrange'
    return variant

def read_cnv_variant(row):
    variant = {}
    variant['HGNC_Symbol'] = row['HGNC_Symbol']
    variant['gene'] = row['HGNC_Symbol']
    variant['Variant'] = row['Variant']
    variant['reported_variant'] = variant['HGNC_Symbol'] + ' ' +variant['Variant']
    if variant['Variant'] == 'Copy Number Gain':
        variant['Variant'] = variant['HGNC_Symbol'] + ' amp'
        variant['pdot'] = 'amp'
    else:
        variant['Variant'] = variant['HGNC_Symbol'] + ' loss'
        variant['pdot'] = 'loss'
    variant['type'] = 'cnv'
    variant['report_status'] = 'not_reported'
    return variant


def read_snv_variant(row):
    variant = {}
    variant['HGNC_Symbol'] = row['HGNC_Symbol']
    variant['Variant'] = row['Variant']
    v = variant['Variant']
    lp = v.find('(')
    if lp != -1:
        rp = v.find(')')
        v = v[lp+1:rp]

    variant['reported_variant'] = variant['HGNC_Symbol'] + ' ' + v
    variant['Quality'] = row['Quality']
    variant['VAF'] = 100.0 * float(row['VAF'])
    chrom = row['Chromosome']
    if isinstance(chrom, numbers.Number):
        chrom = int(chrom)

    variant['chrom'] = str(chrom)
    variant['ref'] = row['Ref']
    variant['alt'] = row['Alt']
    # old pos
    variant['pos_hg19'] = row['Position']
    if not isinstance(variant['pos_hg19'], numbers.Number):
        variant['pos_hg19'] = variant['pos_hg19'].replace(u'\ufeff', '')
    variant['pos_hg19'] = int(variant['pos_hg19'])
    variant['type'] = 'snv'
    variant['report_status'] = 'not_reported'
    return variant



def read_all_variants(patients,path):
    with open(path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            order_id = row['OrderID_External']
            patient = patients[order_id]
            if not row['HGNC_Symbol']== 'NO VARIANTS REPORTED':
                if 'Copy Number' in row['Variant']:
                    cnv = read_cnv_variant(row)
                    patient['unannotated_cnv'].append(cnv)
                elif 'Fusion' in row['Variant']:
                    fusion = read_fusion(row)
                    patient['unannotated_fusion'].append(fusion)
                else:
                    snv = read_snv_variant(row)
                    patient['unannotated_snv'].append(snv)


def read_immune_results_file(path):
    patients = {}
    with open(path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            patient = {'unannotated_snv':[], 'unannotated_cnv':[],'unannotated_fusion':[], 'io_data':{}}
            seen_TMB_Interpretation = False
            for key in row:
                if seen_TMB_Interpretation:
                    patient['io_data'][key] = row[key]
                else:
                    patient[key] = row[key]
                    seen_TMB_Interpretation = 'TMB Interpretation' in key

            patients[patient['OrderID_External']] = patient
    return patients


def n_digit_random(n):
    digits = ''.join(["%s" % randint(0, 9) for num in range(0, n)])
    return digits


def get_fake_order_id(patient,fake_id_dict):
    if patient['OrderID_External'] in fake_id_dict:
        fake_order_id = fake_id_dict[patient['OrderID_External']]
    else:
        fake_order_id = 'P-19-0' + str(n_digit_random(4))
    patient['fake_order_id'] = fake_order_id



def add_patient_data(patient, fake_id_dict,disease_icd_dict,omni_to_jax_disease_dict,disease_path_to_reportable_disease_dict):
    source_dict = {'Lung Cancer': '  lung, biopsy', 'Melanoma': 'skin, biospy','Ovarian Cancer':'biopsy','Uterine Cancer':'biopsy','Breast Cancer':'biopsy','Colorectal Cancer':'biopsy','Unknown Primary Cancer':'biopsy' }
    faker = Faker()

    get_fake_order_id(patient,fake_id_dict)

    genders = ['Male', 'Female']
    disease_path = patient['DiseasePath'].lower()
    if 'ovary' in disease_path or 'endometrial' in disease_path or 'breast' in disease_path:
        patient['sex'] = 'Female'
    else:
        patient['sex'] = faker.word(ext_word_list=genders)

    if patient['sex'] == 'Male':
        patient['name'] = faker.name_male()
    else:
        patient['name'] = faker.name_female()

    dob = faker.date_of_birth(minimum_age=40, maximum_age=90)
    patient['DOB'] = dob.strftime('%m/%d/%Y')

    patient['test_id'] = patient['fake_order_id'] + '-' + str(n_digit_random(3))
    patient['sign_out_date'] = '12/13/2019 10:39 AM ET'
    patient['provider'] = faker.name_male()
    diagnosis = get_ICD(disease_icd_dict, patient)
    if patient['OmniDisease'] in source_dict:
        source = source_dict[patient['OmniDisease']]
    else:
        source = 'biopsy'
    patient['diagnosis'] = diagnosis
    patient['source'] = source
    if patient['DiseasePath'] in disease_path_to_reportable_disease_dict:
        patient['reportable_disease_name'] = disease_path_to_reportable_disease_dict[patient['DiseasePath']]
    else:
        patient['reportable_disease_name'] = patient['OmniDisease']
    add_jax_disease(patient,omni_to_jax_disease_dict)


def add_jax_disease(patient,omni_to_jax_disease_dict):
    if patient['DiseasePath'] in omni_to_jax_disease_dict:
        patient['ckb_diseases'] =[]
        patient['ckb_disease_ids'] = []
        entry = omni_to_jax_disease_dict[patient['DiseasePath']]
        for d in entry:
            patient['ckb_diseases'].append(d['jax_disease_name'])
            patient['ckb_disease_ids'].append(d['jax_disease_id'])
    else:
        patient['ckb_diseases'] = [patient['OmniDisease'].lower()]
        patient['ckb_disease_ids'] = None

def get_ICD(disease_icd_dict, patient):
    if patient['OrderID_External'] in disease_icd_dict:
        diagnosis = disease_icd_dict[patient['OrderID_External']]
    else:
        diagnosis = ''
    return diagnosis


def get_io_base_variant(patient):
    variant = {}
    variant['is_oncogene'] = False
    variant['is_activating'] = False
    variant['is_inactivating'] = False
    variant['is_io'] = True
    variant['type'] = 'io'
    variant['reasons'] = []
    variant['category'] = None
    variant['report_status'] = 'reportable'
    patient['io'].append(variant)
    return variant

def get_pdl1_cutoff_for_disease(patient):
    return 1


def add_io_data(patient):
    ihc_ = patient['PD-L1 IHC']
    if ihc_ == '<1':
        pd_l1_result = 0
    else:
        pd_l1_result = int(ihc_)

    if patient['PD-L1 Clone'] == '22C3':
        patient['pd-l1_primary_result'] = str(patient['PD-L1 IHC']) + '% Tumor Proportion Score (TPS)'
    else:
        patient['pd-l1_primary_result'] = str(patient['PD-L1 IHC']) + '% IC'

    patient['show_pd-l1_secondary'] = patient['PD-L1 Secondary IHC'] != '-'
    if patient['show_pd-l1_secondary']:
        if patient['PD-L1 Secondary Clone'] == '22C3':
            patient['pd-l1_secondary_result'] = str(patient['PD-L1 Secondary IHC']) + '% Tumor Proportion Score (TPS)'
        else:
            patient['pd-l1_secondary_result'] = str(patient['PD-L1 Secondary IHC']) + '% IC'


    variant = get_io_base_variant(patient)
    variant['gene'] = 'PD-L1'
    if pd_l1_result >= get_pdl1_cutoff_for_disease(patient):
        variant['pdot'] = 'positive'
        patient['PD-L1 Interpretation'] = 'positive'
    else:
        variant['pdot'] = 'negative'
        patient['PD-L1 Interpretation'] = 'negative'
    variant['gene_description'] = 'Programmed death-ligand 1 (PD-L1) is a 40kDa type 1 transmembrane protein that plays a major role in suppressing the adaptive arm of immune system. The binding of PD-L1 to the inhibitory checkpoint molecule PD-1 transmits an inhibitory signal based on interaction with phosphatases (SHP-1 or SHP-2) via Immunoreceptor Tyrosine-Based Switch Motif (ITSM) motif. This reduces the proliferation of antigen-specific T-cells in lymph nodes, while simultaneously reducing apoptosis in regulatory T cells (anti-inflammatory, suppressive T cells) - further mediated by a lower regulation of the gene Bcl-2.'

    tmb_interpretation = patient['TMB Interpretation']
    variant = get_io_base_variant(patient)
    variant['gene'] = 'TMB'
    if tmb_interpretation == 'High':
        variant['pdot'] = 'high'
    elif tmb_interpretation == 'Low':
        variant['pdot'] = 'low'
    else:
        variant['pdot'] = 'intermediate'
    variant['gene_description'] = 'Tumor mutational burden is a measurement of mutations carried by tumor cells ' \
                                  'and is a predictive biomarker being studied to evaluate its association with ' \
                                  'response to Immuno-Oncology (I-O) therapy.  Tumor cells with high TMB may have more neoantigens, with an associated increase in cancer-fighting T cells in the tumor microenvironment and periphery. These neoantigens can be recognized by T cells, inciting an anti-tumor response.'


    msi_result = patient['MSI']
    variant = get_io_base_variant(patient)
    variant[
        'gene_description'] = 'Microsatellite instability (MSI) is the condition of genetic hypermutability (predisposition to mutation) that results from impaired DNA mismatch repair (MMR). The presence of MSI represents phenotypic evidence that MMR is not functioning normally.'
    variant['gene'] = 'MSI'
    if msi_result == 'High' or 'Unstable' in msi_result or 'MSI-H' in msi_result:
        variant['pdot'] = 'high'
    elif msi_result == 'Low' or msi_result == 'Stable':
        variant['pdot'] = 'stable'
    else:
        variant['report_status'] = 'not_reported'


def read_summary_interprations(patients,path):
    summary_dict = {}
    with open(path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            summary_dict[row['Order_ID']] = {'likelihood': row['Likelihood'],'molecularsummary': row[
                'MolecularSummary']}
    for order_id in patients.keys():
        patient = patients[order_id]
        patient['likelihood'] = summary_dict[order_id]['likelihood']
        patient['molecularsummary'] = summary_dict[order_id]['molecularsummary']


def read_fake_ids():
    fake_id_dict = {}
    with open('input/name_mapping.csv') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',')
        for row in csv_reader:
            fake_id_dict[row[1]] = row[0]
    return fake_id_dict


def read_disease_icd_dict():
    disease_icd_dict = {}
    with open('input/ICDcodes_for_orders_2020_01_08.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            disease_icd_dict[row['external_order_id']] = row['code'] + ', ' + row['description']
    return disease_icd_dict

def read_omni_to_jax_disease_dict():
    omni_to_jax_disease_dict = {}
    with open('data/tblOS_GLOBAL_JAX_DL_OmniMap.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            disease_entry = {'disease_path': row['DiseasePath'],
                             'jax_disease_id': row['ResourceDiseaseID'],
                             'jax_disease_name': row['ResourceDiseaseName'],
                             'omni_disease': row['OmniDisease'], 'm_code': row['MCode']}
            if row['DiseasePath'] in omni_to_jax_disease_dict:
                entry = omni_to_jax_disease_dict[row['DiseasePath']]
            else:
                entry = []
                omni_to_jax_disease_dict[row['DiseasePath']] = entry
            entry.append(disease_entry)
    return omni_to_jax_disease_dict

def read_omni_reportable_disease_name_dict():
    disease_path_to_reportable_disease_dict = {}
    with open('data/tblPROD_OA_Ref_OmniMap.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            disease_path = row['DiseasePath']
            reportable_name = row['ReportDiseaseName']
            if len(disease_path) > 0 and len(reportable_name)>0:
                disease_path_to_reportable_disease_dict[disease_path] = reportable_name
    return disease_path_to_reportable_disease_dict

def read_variant_groups_dict():
    variant_groups_dict = {}
    with open('data/tblPROD_OA_Ref_VariantGroups.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            vg_entry = {'report_disease_name': row['ReportDiseaseName'],
                             'variant_group': row['VariantGroup'],
                             'variant_name': row['VariantName']}
            if row['ReportDiseaseName'] in variant_groups_dict:
                entry = variant_groups_dict[row['ReportDiseaseName']]
            else:
                entry = []
                variant_groups_dict[row['ReportDiseaseName']] = entry
            entry.append(vg_entry)
    return variant_groups_dict
