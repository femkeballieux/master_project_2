import alminer
import pandas
from astropy.io import ascii
import time

print("Running code")
start = time.time()

# myquery = alminer.target(["AB Aur"])
# alminer.keysearch({"proposal_abstract": ["'high-mass star formation' outflow"]})
# myquery = alminer.keysearch({'target_name': ['GRB021004','SPT0319-47', 'G345']})
# myquery = alminer.keysearch({'proposal_abstract': ['"high-mass star formation" outflow',
#                                                    '"massive star formation" outflow']},
#                             print_targets=False)

observations = alminer.keysearch({'proposal_abstract': ['"high-mass star formation" outflow',
                                    '"massive star formation" outflow']}, print_targets=False)
# print(alminer.explore(observations, allcols=True, allrows=False))
print(alminer.get_info('ang_res_arcsec'))


end = time.time()

print('Takes ', end-start, ' seconds to run')