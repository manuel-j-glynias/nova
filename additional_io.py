import csv
import KMS

lorem = {'description':'Auctor augue mauris augue neque gravida in fermentum et sollicitudin ac orci phasellus egestas tellus rutrum tellus pellentesque eu tincidunt tortor aliquam nulla facilisi cras fermentum odio eu feugiat pretium nibh ipsum consequat in', 'pmid': '8276250'}

io_descriptions = {
'LAG3' :	{'description':'Sustained LAG3 (Lymphocyte activation gene-3) (CD223) expression contributes to a state of lymphocyte exhaustion manifest in impaired proliferation and cytokine production, and shows a remarkable synergy with PD-1 to inhibit immune responses.', 'pmid': '28258692'},
'TIM3' :	{'description':'TIM-3 (T-cell immunoglobulin mucin-3, also known as HAVCR2)  inhibits antitumor immunity by mediating T-cell exhaustion.  Blocking the TIM-3 pathway enhances cancer immunity and increases the production of interferon-gamma (IFN-γ) in T cells.', 'pmid':'30410357'},
'VISTA' :	{'description':'VISTA (V-domain Ig suppressor of T cell activation) is a negative immune checkpoint protein, which suppress T cell activation.VISTA overexpression on tumor cells interferes with protective antitumor immunity in the host which is dependent on T cells.', 'pmid':'28258699'},
'CD137' :	{'description':'CD137 (4-1BB/TNFRSF9)  is a co-stimulatory molecule belonging to the tumor necrosis factor receptor superfamily.  Agonistic anti-CD137 mAb has been demonstrated to induce and improve anti-cancer immunity in several models of cancer.', 'pmid':'31013788'},
'CD27' :	{'description':'CD27 (TNFRSF7)is a co-stimulatory molecule belonging to the tumor necrosis factor receptor superfamily. CD27 agonists (either CD70 or agonist anti-CD27 mAb) costimulate human T cells to undergo proliferation and cytokine production in vitro.  CD27-driven CD8+ T-cell activation can occur directly without a requirement for CD4+ T cells.', 'pmid':'291180060'},
'CD40' :	{'description':'Cluster of differentiation 40 (CD40) is a cell surface molecule of the tumour necrosis factor receptor family. CD40 is expressed on antigen-presenting cells (APCs), for example, dendritic cells and myeloid cells, and is critical for their activation and proliferation.  Administration of an agonistic antibody directed against CD40 produced protective T cell immunity in murine models of cancer.', 'pmid':'27927088'},
'GITR' :	{'description':'Glucocorticoid-induced tumor necrosis factor receptor-related protein (GITR) promotes effector T cell  functions and hamper regulatory T cell suppression. Although GITR agonism is not sufficient to activate cytolytic T cells due to persistent exhaustion, T cell reinvigoration with PD-1 blockade can overcome resistance of  tumors to anti-GITR monotherapy.', 'pmid':'31036879'},
'ICOS' :	{'description':'The Inducible T cell Costimulator (ICOS, CD278) is a receptor in the CD28 family of B7-binding proteins expressed mostly by activated T cells. Upon binding to its ligand ICOS-L expressed on Antigen Presenting Cells, T cells are further activated by ICOS, resulting in their expansion and production of effector cytokines.', 'pmid':'29468927'},
'OX40' :	{'description':'OXO40 is a co-stimulatory molecule belonging to the tumor necrosis factor receptor superfamily.  In vitro studies have shown that stimulation of OX40 enhances proliferation and expression of effector molecules and cytokines by human T cells.  The effects on CD8+ T cells are largely indirect, being attributable to enhanced T-cell help provided by the CD4+ population.', 'pmid':'291180060'},
'CCL2' :	{'description':'CCL2, also named monocyte chemotactic protein-1, is a potent chemokine for the recruitment of tumor-associated macrophages (TAMs)  via activating its receptor CCR2 in the process of tumor invasion and metastasis. Growing preclinical studies showed that targeting TAMs by blocking the CCL2/CCR2 axis inhibited the tumor development and might be a potential therapeutic revenue for several malignancies.', 'pmid':'31024838'},
'CCR2' :	{'description':'CCL2, also named monocyte chemotactic protein-1, is a potent chemokine for the recruitment of tumor-associated macrophages (TAMs)  via activating its receptor CCR2 in the process of tumor invasion and metastasis. Growing preclinical studies showed that targeting TAMs by blocking the CCL2/CCR2 axis inhibited the tumor development and might be a potential therapeutic revenue for several malignancies.', 'pmid':'31024838'},
'CSF1R' :	{'description':'Colony stimulating factor-1 (CSF1R) receptor, a memeber of the type III growth factor receptor family, regulates the differentiation and survival of macrophages. Blockade of CSF1 or its receptor has shown to markedly decrease the infiltration of macrophages at the tumor site, and to attenuate primary tumor growth, reduce metastatic potential and improve long-term survival of some tumor-bearing mice.', 'pmid':'31215373'},
'ADORA2A' :	{'description':'The adenosine A2A receptor, also known as ADORA2A, is an adenosine receptor, has also been shown to play a regulatory role in the adaptive immune system.  Tumor cells and tumor-associated Treg adenosine generation to dampen intratumoral immune responses, through the suppression of T cell effector functions via stimulation of the ADORA2A on effector T cells .', 'pmid':'31024543'},
'CD39' :	{'description':'CD39 is an ecto-nucleoside triphosphate diphosphohydrolase which hydrolyzes extracellular ATP to generate generate immunosuppressive adenosine nucelosides. Targeting immunosuppressive adenosine by inhibiting CD39 may restore anti-tumor responses or boost the efficacy of other anti-cancer therapies.', 'pmid':'28258700'},
'IDO1' :	{'description':'Indoleamine 2,3-dioxygenase 1 (IDO1) is a principle enzyme in Trp catabolism, and one mediator of immunosuppression in T cell–inflamed tumors is the tryptophan–kynurenine–aryl hydrocarbon receptor (Trp–Kyn–AhR) pathway.  Elevated activity of this pathway is associated with increased tumor grade and poor prognosis in many cancers.  Gene expression of IDO1 is strongly correlated with the expression of PD-1.', 'pmid':'30377198'},
'TGFB1' :	{'description':'TGFB1 ( transforming growth factor β 1) attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Combination TGFβ inhibition and immunotherapy induces complete responses in mouse models.', 'pmid':'29443960'},
}

#https://www.ncbi.nlm.nih.gov/pubmed/29443960/
io_recommendations = {
    'GITR': 'Glucocorticoid-induced TNFR-related protein (GITR) is a surface protein that is upregulated upon T cell '
            'activation.  In in vivo models, administration of a GITR agonist antibody is associated with reduction of intratumoral Treg accumulation and potentiation of antitumor CD8+ effector T cell function (Cohen, Schaer, Liu, et al., 2010; Ko et al., 2005; Mitsui et al., 2010), as well as antitumor activity (Cohen et al., 2010; Ko et al., 2005; Turk et al., 2004). When given in combination with PD-1 blockade, increased activity was also seen. ',
    'CD137' : 'Preclinical studies in mice demonstrated strong antitumor response with agonist anti-CD137 mAbs ('
              'Melero et al., 1997). In fact, stimulation of CD137 is one of the most powerful antibody-based cancer immunotherapeutic strategies in mouse models (Vinay and Kwon, 2012; Wilcox et al., 2002). In addition, agonist anti-CD137 synergizes with radiotherapy and chemotherapy (Ju et al., 2008; Shi and Siemann, 2006).',
    'TIM3' : 'TIM-3, a co-inhibitory molecule, marks exhausted T cells in the tumor microenvironment and during '
             'chronic viral infection, and is often coexpressed with PD-1 on tumor-antigen-specific T cells in cancer patients. Prolonged T-cell exhaustion results in the inability to respond to anti-PD1.  Clinically, TIM-3 is also upregulated in anti-PD1-resistant tumors. TIM-3 inhibition could overcome anti-PD-1/L1 resistance driven by TIM-3 mediated exhaustion and provides the opportunity to enhance antitumor T-cell immunity',
        'CD39': 'Combination of a PD-1 axis inhibitor with a CD39/ADORA2A antagonist. CPI-444 (from Corvus Pharmaceuticals) has been tested in 47 adults with advanced solid tumors of varying histology who had previously failed at least 5 lines of treatment. With a median follow-up time of 12 weeks (range 2-32), the overall disease control rate (CR, PR or SD) was 45% per RECIST. The proportion of patients with disease control were similar for pts treated with CPI-444 alone and for patients treated with CPI-444 plus atezolizumab',
        'ADORA2A': 'Combination of a PD-1 axis inhibitor with a CD39/ADORA2A antagonist. CPI-444 (from Corvus Pharmaceuticals) has been tested in 47 adults with advanced solid tumors of varying histology who had previously failed at least 5 lines of treatment. With a median follow-up time of 12 weeks (range 2-32), the overall disease control rate (CR, PR or SD) was 45% per RECIST. The proportion of patients with disease control were similar for pts treated with CPI-444 alone and for patients treated with CPI-444 plus atezolizumab.',
        'IDO1': 'Recent results from the Phase 3 ECHO-301/KEYNOTE-252 study combining the IDO1 inhibitor epacadostat with pembrolizumab in advanced melanoma did not show improved survival. This compares to phase 2 results for the IDO1 inhibitor BMS-986205 in combination with nivolumab presented at AACR in 2017 that showed ORR of 20% to 40% in heavily pretreated cervical and bladder cancer.',
        'LAG3' : 'Encouraging clinical results were reported at ESMO 2017 for the combination of nivolumab and a LAG3 inhibitor in patients with metastatic melanoma. Of 68 heavily pretreated patients, including prior immunotherapy, there was 1 complete response and 6 partial responders, with no association to PD-L1 expression but an excellent correlation to LAG3 expression.',
        'ICOS' : 'At ASCO 2018 phase 1 & 2 results were reported for the ICOS agonist JTX-2011 with and without combination nivolumab. In 75 patients in the phase 2 portion combination of this trial the disease control rate (CR, PR, SD) was 32%. Of note, increased CD4+ ICOS high T cells were observed in 7/7 (n=29) evaluable patients with PR.',
        'CD27' :'For the CD27 agonist varililumab a single phase 1a single agent study reported 1 partial response and 8 stable disease in 25 patients (PMID: 28463630).',
        'CD40': 'In 2015 an abstract was reported at AACR for 24 melanoma patients treated with a combination of a CD40 agonist and CLTA-4 inhibitor. The overall clinical benefit was not beyond what was expected for single agent CTLA-4 inhibitor.',
        'OX40' : 'A single agent phase 1 trial (PMID: 24177180) of an anti-OX40 monoclonal antibody showed no complete or partial responses in 30 patients of multiple histology, although stable disease was observed in 17 patients. A second phase1 human clinical trial of another anti-OX40 monoclonal antibody, BMS-986178, in combination with nivolumab of multiple tumor types was presented at SITC 2017, but response rates were not disclosed.',
        'VISTA' : 'For the dual PD-1 and VISTA inhibitor CA-170, at SITC 2017 the first in human phase 1a results were presented in 39 patients with advanced cancer of multiple histology of which the majority were checkpoint inhibition naive. The results were not overly encouraging with best response of stable disease in 8 patients.',
        'CSF1R' : 'Among the various CSF1R or CSF1 antagonists only five, including PLX3397, JNJ-40346527, cabiralizumab,  emactuzumab, and AMG820 have reported only phase 1 toxicity studies with minimal to modest adverse events. It should be noted that the majority of these studies were completed in either tenosynovial giant cell tumors or classical Hodgkin lymphoma, with only AMG820 tested in advanced solid tumors. In this study best response was SD in 8 pts (32%) (NCT01444404).',
        'TGFB1' : 'For the TGFB1 antagonist galunisertib, results have been published as single agent therapy primarily in somewhat uncommon tumor types such as glioma or hepatocellular carcinoma with modest results (ORR 5-10% or less), but these studies to date were not in combination with a PD-1 axis inhibitor which would likely improve the results.',
        'CCL2' : 'Combination of chemotherapy, a checkpoint inhibitor, and CCR2 anatagonist (BMS-813160). There have '
                 'been encouraging clinical efficacy results for CCR2 inhibition reported in pancreatic cancer in '
                 'combination with FOLFIRINOX chemotherapy (PMID: 27055731). In this study 47 patients with '
                 'pancreatic ductal adenocarcinoma were treated with either FOLFIRINOX alone (n=8) or with FOLFIRINOX plus PF-04136309 (n=39), CCR2 antagonist. In the combination arm 16 (49%) of 33 patients receiving FOLFIRINOX plus PF-04136309 who had undergone repeat imaging achieved an objective tumour response, compared to none in the FOLFIRINOX alone group.',
        'CCR2' : 'Combination of chemotherapy, a checkpoint inhibitor, and CCR2 anatagonist (BMS-813160). There have '
                 'been encouraging clinical efficacy results for CCR2 inhibition reported in pancreatic cancer in '
                 'combination with FOLFIRINOX chemotherapy (PMID: 27055731). In this study 47 patients with '
                 'pancreatic ductal adenocarcinoma were treated with either FOLFIRINOX alone (n=8) or with FOLFIRINOX plus PF-04136309 (n=39), CCR2 antagonist. In the combination arm 16 (49%) of 33 patients receiving FOLFIRINOX plus PF-04136309 who had undergone repeat imaging achieved an objective tumour response, compared to none in the FOLFIRINOX alone group.'
        }

IO_cut_off = 75
reportable_io_markers = ['CD8','TIM3','CCR2','GITR','VISTA','CSF1R','TNF','CD38','IDO1','TGFB1','CTLA4','LAG3','CD40','ADORA2A','OX40','CD137','CD27','ICOS']

def get_io_description(key):
    value = lorem
    if key in io_descriptions:
        value = io_descriptions[key]
    return value

def get_io_recommendation(key):
    value = lorem['description']
    if key in io_recommendations:
        value = io_recommendations[key]
    return value


io_markers = ['CD137','CD27','CD28','CD40','CD40LG','CD80','CD86','GITR','GZMB','ICOS','ICOSLG','IFNG','OX40','OX-40L','TBX21',
             'CXCL10','CXCR6','DDX58','GATA3','IL10','IL1B','MX1','STAT1','TGFB1','TNF',
              'CD2','CD20','CD3','CD4','CD8','FOXP3','KLRD1','SLAMF4',
              'BTLA','CTLA4','LAG3','PD-1','PD-L1','PD-L2','TIM3','VISTA','TNFRSF14',
              'ADORA2A','CCL2','CCR2','CD163','CD38','CD39','CD68','CSF1R','IDO1'
]

def make_trial(io,trial):
    s = ''
    for drug in trial['drugs']:
        if len(s) > 0:
            s += ' + '
        s += drug
    availability = 'Available Trial(s): ' + trial['nct_id'] + ' (' + trial['phase']
    for loc in trial['locations']:
        availability += ', ' + loc['facility_name']
    availability += ')'

    fake_trail = {'therapy':s, 'availability':availability}
    fake_trail['efficacy'] = get_io_recommendation(io)
    return fake_trail


def get_trial_io(patient,io,drug):
    zip = '14203'
    clincal_trial_distance = 150

    trials = KMS.match_go_local_trials_by_drug(drug,zip,clincal_trial_distance,patient['DOB'])
    return trials



def add_additional_io(patient,io_drug_dict):
    patient['additional_io_details'] = []
    for io in io_markers:
        if io in reportable_io_markers and not patient['io_data'][io]=='-':
            if int(patient['io_data'][io]) >= IO_cut_off or (io=='CD8' and int(patient['io_data']['CD8']) <= 25):
                if io in io_drug_dict:
                    drug_list = io_drug_dict[io]
                    for d in drug_list:
                        trials = get_trial_io(patient,io,d)
                        if len(trials)>0:
                            trial_list =[]
                            for trial in trials:
                                trial_list.append(make_trial(io,trial))
                            detail = {'marker':io, 'rank': patient['io_data'][io], 'interpretation':'High', 'description':get_io_description(io),
                                      'trials': trial_list}
                            patient['additional_io_details'].append(detail)



def get_io_drug_dict():
    file_path = 'data/tblPROD_OA_OA_MarkerReporting.csv'
    io_drug_list = []
    with open(file_path, newline='') as io_drugs:
        io_reader = csv.DictReader(io_drugs)
        for io_drug in io_reader:
            if len(io_drug['Therapies']) > 0:
                io_drug['drug_list'] = io_drug['Therapies'].split(";")

            io_drug_list.append(dict(io_drug))

    io_drug_dict = {}
    for io_drug in io_drug_list:
        if 'drug_list' in io_drug:
            io_drug_dict[io_drug['Marker']] = io_drug['drug_list']
    return io_drug_dict

