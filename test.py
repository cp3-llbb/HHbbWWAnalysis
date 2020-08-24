file_ = 'jobs.log'
with open(file_) as fp:
    logfiles = fp.readlines()
    for lfile_ in logfiles:
        print(lfile_.strip())
        with open(lfile_.strip()) as log:
            loglines = log.readlines()
            print([item.strip() for item in loglines if 'Batch node' in item])
 
