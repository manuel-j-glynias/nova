import numbers
import csv
from random import randint
from faker import Faker

def read_fusion(row):
    variant = {}
    variant['HGNC_Symbol'] = row['HGNC_Symbol']
    variant['gene'] = row['HGNC_Symbol']
    variant['Variant'] = row['Variant']
    variant['pdot'] = 'rearrange'


def read_cnv_variant(row):
    variant = {}
    variant['HGNC_Symbol'] = row['HGNC_Symbol']
    variant['gene'] = row['HGNC_Symbol']
    variant['Variant'] = row['Variant']
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


def add_patient_data(patient):
    disease_dict = {'Lung Cancer': ' C34, Malignant neoplasm of bronchus and lung',
                    'Melanoma': 'C43.59, Malignant melanoma of other part of trunk, Unknown',
                    'Ovarian Cancer':'C56.9, Malignant neoplasm of unspecified ovary',
                    'Uterine Cancer':'C54.1, Malignant neoplasm of endometrium',
                    'Breast Cancer':'C50.919, Malignant neoplasm of unspecified site of unspecified female breast',
                    'Colorectal Cancer':'C18.9, Malignant neoplasm of colon, unspecified',
                    'Unknown Primary Cancer':'C80.1, Malignant (primary) neoplasm, unspecified'
                    }
    source_dict = {'Lung Cancer': '  lung, biopsy', 'Melanoma': 'skin, biospy','Ovarian Cancer':'biopsy','Uterine Cancer':'biopsy','Breast Cancer':'biopsy','Colorectal Cancer':'biopsy','Unknown Primary Cancer':'biopsy' }
    faker = Faker()


    fake_order_id = 'P-19-0' + str(n_digit_random(4))
    test_id = fake_order_id + '-' + str(n_digit_random(3))
    patient['fake_order_id'] = fake_order_id
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

    patient['test_id'] = test_id
    patient['sign_out_date'] = '12/13/2019 10:39 AM ET'
    patient['provider'] = faker.name_male()
    if patient['OmniDisease'] in disease_dict:
        diagnosis = disease_dict[patient['OmniDisease']]
    else:
        diagnosis = ''
    if patient['OmniDisease'] in source_dict:
        source = source_dict[patient['OmniDisease']]
    else:
        source = 'biopsy'
    patient['diagnosis'] = diagnosis
    patient['source'] = source
    patient['ckb_disease'] = patient['OmniDisease'].lower()


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
    pd_l1_result = int(patient['PD-L1 IHC'])

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