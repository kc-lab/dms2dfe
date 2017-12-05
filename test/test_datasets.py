def test_datasets():
    import subprocess
    import os
    from os.path import exists

    def listfiles(path):
        files = []
        for dirName, subdirList, fileList in os.walk(path):
            dir = dirName.replace(path, '')
            for fname in fileList:
                if not '.png' in fname: 
                    files.append(os.path.join(dir, fname))
        files=[f for f in files if (('/data_fit/' in f) or ('/data_comparison/' in f) or ('/plots/' in f))]
        return files

    repo_n='ms_datasets'
    datasets=[
              # 'Firnberg_et_al_2014',
              'Olson_et_al_2014',
              # 'Melnikov_et_al_2014',
              ]
    datasets_dh='%s' % (repo_n)
    if not exists(datasets_dh):
        subprocess.call('git clone https://github.com/rraadd88/%s.git' % (repo_n),shell=True)

    os.chdir('%s/analysis' % (repo_n))
    for prj_dh in datasets:
        print prj_dh
        if exists(prj_dh):
            if not exists(prj_dh+"/plots"):
                subprocess.call('dms2dfe %s' % (prj_dh),shell=True)
            x = listfiles(prj_dh)
            y = listfiles('../outputs/%s' % prj_dh)
            print(">files in analysis directory:",x)
            print(">files in outputs directory:",y)
            files_only_in_either = set(x) ^ set(y)
            if (len(y)!=0) or (len(files_only_in_either)==0):
                print ">>> SUCCESS"
            else:
                print ">>> TEST NOT SUCCESSFUL. Following files are different betn analysis and output"
                print files_only_in_either
        break
test_datasets()
