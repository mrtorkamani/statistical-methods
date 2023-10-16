import numpy as np
from scipy.special import hyp2f1


class SNLSLikelihoodModule(object):

    def __init__(self, snlsdata_path=None, intrinsic_disp = [0.0675, 0.1133, 0.0815, 0.0989]):
        """
        Class for calculating SNLS likelihood.

        :param snlsdata_path (optional): path to folder containing SNLS data, default: None (i.e. in active folder)
        """
        self.lcdparams = self.load_data(snlsdata_path)
        self.nsn = self.lcdparams.shape[0]
        self.intrinsic_disp = np.array(intrinsic_disp)

    def load_data(self, snlsdata_path):
        """
        Routine for loading the SNLS data

        :param snlsdata_path: path to folder containing SNLS data
        :returns lcdparams: light curve data of SNLS supernovae
        """
        if snlsdata_path == None:
            snlsdata_path = ''
        lcdparams = np.genfromtxt(snlsdata_path+'snls_3rdyear_lcparams.txt', names = True, dtype=None, encoding=None)
        return lcdparams

    def getLCDMLuminosityDistance(self, z, omegam, zhel = None):
        """
        Calculate the Hubble-constant free luminosity distance (i.e. in units of MPc/H0)
        using scipy.special.hyp2f1 for a flat LambdaCDM cosmology.

        :param omegam: Matter density today
        :param z: Redshift or array of redshifts
        :returns dL: luminosity distance to redshifts in z
        """
        if zhel is None:
            zhel = z
        omoverol = omegam / (1 - omegam)
        out = (1+z)*hyp2f1(1./3., .5, 4./3., -(1+z)**3 * omoverol)
        out -= hyp2f1(1./3., .5, 4./3., -omoverol)
        return out / np.sqrt(1 - omegam) * (zhel + 1)

    def getModelMag(self, omegam, alpha, beta, M):
        """
        Routine for calculating model magnitudes of all supernovae in SNLS dataset

        :param omegam: Matter density today
        :param alpha: stretch parameter alpha
        :param beta: color parameter beta
        :param M: absolute magnitude parameter M
        :returns mmod: model magnitudes for SNLS supernovae
        """
        dL = self.getLCDMLuminosityDistance(self.lcdparams['zcmb'], omegam, self.lcdparams['zhel'])
        mmod = 5 * np.log10(dL) - alpha * (self.lcdparams['s'] - 1) + beta * self.lcdparams['color'] + M
        return mmod

    def getVariance(self, alpha, beta):
        """
        Routine for calculating the independent variance of each supernova magnitudes in SNLS data.

        :param alpha: stretch parameter alpha
        :param beta: color parameter beta
        :returns Var: independent variance of all supernovae magnitudes in SNLS data
        """
        Var = (self.lcdparams['dmb']*self.lcdparams['dmb'] + (alpha*alpha) * self.lcdparams['ds'] * self.lcdparams['ds']
               + (beta*beta) * self.lcdparams['dcolor'] * self.lcdparams['dcolor']
               + (2*alpha) * self.lcdparams['cov_m_s'] - (2*beta) * self.lcdparams['cov_m_c']
               - (2*alpha*beta) * self.lcdparams['cov_s_c'] + self.intrinsic_disp[self.lcdparams['set']]
               + (5*5/np.log(10)/np.log(10)) * (self.lcdparams['dz']*self.lcdparams['dz']/self.lcdparams['zcmb']
                                          /self.lcdparams['zcmb']))
        return Var


    def getLogLikelihood(self, omegam, alpha, beta, M):
        """
        Routine for evaluating the log-likelihood.

        :param omegam: Matter density today
        :param alpha: stretch parameter alpha
        :param beta: color parameter beta
        :param M: absolute magnitude parameter M
        """
        if omegam > 1:
            return -1e30
        else:
            mmod = self.getModelMag(omegam, alpha, beta, M)
            Var = self.getVariance(alpha, beta)
            N = np.log(Var).sum() + self.nsn * np.log(2 * np.pi)
            chi2 = ((self.lcdparams['mb'] - mmod) * (self.lcdparams['mb'] - mmod) / Var).sum()
        return -.5 * (chi2 + N)
