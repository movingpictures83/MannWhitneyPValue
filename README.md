# MannWhitneyPValue
# Language: Python
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 2.0, Python 3.6
# Dependency: scipy_1.4.1, pandas_1.1.5

PluMA plugin that runs the Mann-Whitney test (Mann and Whitney, 1947)
and computes P-Values to determine differentiating taxa for two
sets of samples, at two different levels.

This is useful for example when you have samples with mixed properties
(i.e. diseased and healthy, high-fat vs. low-fat diet).

The plugin will find differentiating taxa between one level (i.e. diseased and healthy),
then within each find differentiating taxa between the other level (i.e. high-fat vs. low-fat diet).

The former is 'coarse', the latter is 'fine'.  Output files produced
will be <prefix>.coarse.txt and <prefix>.fine.txt

Input is a tab-delimited TXT file of parameters:
metadata: CSV file of samples and groups
abundance: CSV file of abundances
group1: General group 1
group2: General group 2
group1a: Sub-group 1
group2a: Sub-group 2

