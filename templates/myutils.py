class Mytempclass:
    def __init__(self, x,mean,median,std):
        self.mean=mean
        self.median=median
        self.std=std
        self.x=x

    def tempfuncy(self):
        from scipy.interpolate import interp1d

        tempfuncy = interp1d(self.x,self.mean, kind='cubic', bounds_error=False)
        if np.sum(np.isnan(tempfuncy(self.x)))==len(self.x) or np.std(np.array(tempfuncy(self.x)[~np.isnan(tempfuncy(self.x))]))>10:
            tempfuncy = interp1d(self.x,self.mean, kind='linear', bounds_error=False)
        return tempfuncy

    def tempfuncstd(self):
        from scipy.interpolate import interp1d
        tempfuncstd = interp1d(self.x,self.std, kind='cubic', bounds_error=False)
        if np.sum(np.isnan(tempfuncstd(self.x)))==len(self.x) or np.std(np.array(tempfuncstd(self.x)[~np.isnan(tempfuncstd(self.x))]))>10:
            tempfuncstd = interp1d(self.x,self.std, kind='linear', bounds_error=False)
        return tempfuncstd


