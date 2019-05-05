import math
import pandas as pd

from statsmodels.sandbox.stats.multicomp import multipletests

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from gdcapi.endpoint.cases import Cases

def get_survival(cancer_type='TCGA-LIHC'):
    cases_endpoint = Cases()
    cases = cases_endpoint(cases_endpoint.project.project_id == cancer_type,
                           fields = [ cases_endpoint.submitter_id,
                                      cases_endpoint.diagnoses.days_to_death,
                                      cases_endpoint.diagnoses.days_to_last_follow_up,
                                      cases_endpoint.diagnoses.vital_status],
                           size=1000)  # TODO handle size cancer type dependent, 1000 should work for all though
    vital_status = [vs[0] if vs else 'NA' for vs in cases.diagnoses.vital_status()]
    days_to_death = [dtd[0] if dtd else math.nan for dtd in cases.diagnoses.days_to_death()]
    days_to_last_followup = [dtlfup[0] if dtlfup else math.nan
                             for dtlfup in cases.diagnoses.days_to_last_follow_up()]
    participants = [p for p in cases.submitter_id()]
    survival = pd.DataFrame({
                   'participant': participants,
                   'vital_status': vital_status,
                   'days_to_death': days_to_death,
                   'days_to_last_followup': days_to_last_followup})
    if set(survival.vital_status.unique()) != {'dead', 'alive'}:
        pass
        # TODO do something, for LIHC its fine
    survival['dead'] = [1 if survival.iloc[i, 3] == 'dead' else 0 for i in range(len(survival))]
    survival['days_to_event'] = [survival.iloc[i, 0] if not math.isnan(survival.iloc[i, 0])
                                                     else survival.iloc[i, 1]
                                 for i in range(len(survival))]
    survival.dropna(subset=['days_to_event'], inplace=True)
    return survival

def simple_survival_analysis(participants, survival, label, full_submitter_id=True):
    if not full_submitter_id:
        participants = [ind for ind in survival['participant'] if ind[-4:] in participants]
    T = survival['days_to_event']
    E = survival['dead']
    ix = survival['participant'].isin(participants)
    kmf = KaplanMeierFitter()
    kmf.fit(T[~ix], E[~ix], label='others')
    ax = kmf.plot()
    kmf.fit(T[ix], E[ix], label=label)
    kmf.plot(ax=ax)
    results = logrank_test(T[ix], T[~ix], event_observed_A=E[ix], event_observed_B=E[~ix])
    return results

def two_group_survival_analysis(survival, participants1, participants2, label1, label2, full_submitter_id=False):
    if not full_submitter_id:
        participants1 = [ind
                         for ind in survival['participant'] if ind[-4:] in participants1]
        participants2 = [ind for ind in survival['participant'] if ind[-4:] in participants2]
    T = survival['days_to_event']
    E = survival['dead']
    ix1 = survival['participant'].isin(participants1)
    ix2 = survival['participant'].isin(participants2)
    kmf = KaplanMeierFitter()
    kmf.fit(T[ix1], E[ix1], label=label1)
    ax = kmf.plot()
    kmf.fit(T[ix2], E[ix2], label=label2)
    kmf.plot(ax=ax)
    results = logrank_test(T[ix1], T[ix2], event_observed_A=E[ix1], event_observed_B=E[ix2])
    return results


def survival_analysis(participants, survival, label, full_submitter_id=True):
    if not full_submitter_id:
        participants = [ind for ind in survival['participant'] if ind[-4:] in participants]
    T = survival['days_to_event']
    E = survival['dead']
    ix = survival['participant'].isin(participants)
    kmf = KaplanMeierFitter()
    kmf.fit(T[~ix], E[~ix], label='others')
    #ax = kmf.plot()
    kmf.fit(T[ix], E[ix], label=label)
    #kmf.plot(ax=ax)
    results = logrank_test(T[ix], T[~ix], event_observed_A=E[ix], event_observed_B=E[~ix])
    return results
    
def survival_all_genewise(subgraphs, survival, threshold=30):
    node_histogram = subgraphs.node_histogram()
    node_histogram = node_histogram[node_histogram['count'] > threshold]
    results = {}
    for gene in node_histogram['gene']:
        subgs = subgraphs.with_any_of(gene)
        results[gene] = survival_analysis(subgs, survival, '')
    genes = list(results.keys())
    corrected_pvals = list(multipletests([results[gene].p_value for gene in genes], method='fdr_bh')[1])
    corrected_pvals = dict(zip(genes, corrected_pvals))
    return results, corrected_pvals

def get_survival_y(survival):
    data = survival[['dead', 'days_to_event', 'participant']].copy()
    data.set_index('participant', inplace=True)
    y = zip(data['dead'], data['days_to_event'])
    y = [(bool(not_censored), int(time)) for not_censored, time in y]
    participants = list(data.index)
    return y, participants