import gefes

#creation of a binner

oil = gefes.projects['oilsand']
oil.binner.new("concoct_5000",{'type' : 'GefesCONCOCT', 'args' : {'max_clusters' : 100, 'max_freq' : 0.15, 'min_length' : 5000, 'transform' : None}} )


#running it once

oil.binner['concoct_5000'].cluster()

#if binner not loaded, load it
oil.binner['concoct_5000'].load()

#run  bin pipeline
for b in oil.binner['concoct_5000']: b.runner.run_slurm()


