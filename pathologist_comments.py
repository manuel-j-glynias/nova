from jinja2 import Environment, BaseLoader
from datetime import datetime
import re



def get_age_from_birthdate(patient):
    dob = patient['DOB']
    datetime_object = datetime.strptime(dob, '%m/%d/%Y')
    new_date = datetime.today() - datetime_object
    age = int(new_date.days / 365.25)
    return age

def num_strongly_recommended(patient):
    num = 0
    # recs = patient['recommendations']
    # if 'strongly_recommended' in recs:
    #     num = len(recs['strongly_recommended'])
    return num

def any_caution(patient):
    num = 0
    # recs = patient['recommendations']
    # if 'caution' in recs:
    #     num = len(recs['caution'])

    return num > 0


# def get_catution(patient):
#     dict = patient['recommendations']['caution']
#     keys = list(dict.keys())
#     key = keys[0]
#     entry = dict[key][0]
#     variant = entry['variant']
#     drug = entry['therapy']
#     return 'Although the patient is ' + variant + ', ' + entry['efficacy_evidence'];

def starts_with_vowel(word):
    vowel = word[0].lower() in "aeiou"
    return vowel


def get_class_and_drug_phrase(patient, strength,index):
    pretty_drug_class = {'egfr inhibitor 3rd gen':'third generation EGFR inhibitor',
                         'pd-l1/pd-1 antibody':'immune checkpoint inhibitor',
                         'ctla4 antibody': 'CTLA-4 antibody'}
    dict = patient['recommendations'][strength]
    keys = list(dict.keys())
    key = keys[index]
    entry = dict[key][0]
    drug_class = entry['drug_classes']
    drug = entry['therapy']
    dc = drug_class[0]
    if dc in pretty_drug_class:
        dc = pretty_drug_class[dc]
    if '+' in drug:
        dc2 = drug_class[1]
        if dc2 in pretty_drug_class:
            dc2 = pretty_drug_class[dc2]
        dc += '/' + dc2 + ' combination'
    article = 'A'
    if index > 0:
        article = 'a'
    if starts_with_vowel(dc):
        article += 'n'
    phrase = article + ' ' + dc + ' such as ' + drug
    return phrase

def get_interpretation():
    interpretation = ''
    with open('templates/interpretation.html', 'r') as file:
        interpretation = file.read().replace('\n', '')
        interpretation = re.sub('\s+', ' ', interpretation).strip()
    return interpretation


def execute_interpretation(patient, interpretation):
    output = ''
    if interpretation is not None:
        template = Environment(loader=BaseLoader).from_string(interpretation)
        template.globals['get_age_from_birthdate'] = get_age_from_birthdate
        template.globals['num_strongly_recommended'] = num_strongly_recommended
        template.globals['any_caution'] = any_caution
        template.globals['get_class_and_drug_phrase'] = get_class_and_drug_phrase
        # template.globals['get_catution'] = get_catution

        output = template.render(patient=patient)
    return output

def get_draft_comments(patient):
    interpretation = get_interpretation()
    output = execute_interpretation(patient, interpretation)
    return output


