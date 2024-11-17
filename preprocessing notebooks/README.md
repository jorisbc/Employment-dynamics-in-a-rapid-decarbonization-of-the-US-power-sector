Execute the files in the following order:

0) Supply Use tables to symmetric IO tables calculation.ipynb
- calculate the (domestic) supply and use tables
1) Crosswalk BLS to OMN
- Construction of crosswalk
2) Read and load BLS industry-occupation data
- Reads relevant industries and occupations from BLS
3) Construct supra-adjacency matrix: use full tables (not only domestic)
- Combines data and crosswalk
4) crosswalk ASEC to BLS and vice versa.ipynb
- crosswalk for occupational mobility network
5) Utilities workers disaggregated
- Includes disaggregated electricity sector workers in data
6) ASEC occ mob network crosswalk.ipynb
- transform the occupational mobility network from ASEC codes to BLS occupation codes
