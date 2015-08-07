
def main():

    import numpy as np

    data = np.loadtxt('/usr/users/ezbc/research/scripts/multicloud/analysis/clouds/cloud_errors.csv',
                      skiprows=1, delimiter=',')

    rms = np.mean(data[:,0]**2)**0.5
    std = np.std(data[:,0])

    print rms, std


if __name__ == '__main__':
    main()

