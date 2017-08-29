InputFile = './combined_wCut.3.txt'
InputFile2 = './combined_wCut.3_upweight.txt'
OutputFile = 'galaxy_bias_combined_wCut.3.pdf'
import numpy as np
from matplotlib import pyplot as plt
z = .85
om = .31
xi0 = np.loadtxt('./xi0isoChallenge_matterpower6.0.dat').transpose()[0]
xdata = np.loadtxt(InputFile).transpose()[0]
def xigg(x,b):
    import sys
    sys.path.append('/global/homes/h/huikong/eboss/LSSanalysis')
    from Cosmo import distance
    d = distance(om,1.-om)
    D = d.D(z)
    f = d.omz(z)**.557
    xi = np.loadtxt('./xi0isoChallenge_matterpower6.0.dat').transpose()[1]
    xi = xi/(1.+2/3.*.4+0.2*.4**2.)
    xi = b**2.*(1.+2/3.*f/b + 0.2*(f/b)**2)*xi*D
    y = np.zeros(len(x))
    for i in range(0,len(x)):
        for j in range(0,len(xi0)):
           if x[i] <= xi0[j]:
              if j == 0:
                  y[i] = xi[j]*x[i]/xi0[j]
                  y[i] = y[i]*x[i]*x[i]
                  break
              else:
                  y[i] = xi[j-1]*(xi0[j]-x[i])/(xi0[j]-xi0[j-1])+xi[j]*(x[i]-xi0[j-1])/(xi0[j]-xi0[j-1])
                  y[i] = y[i]*x[i]*x[i]
                  break                
    return y

def xigg2(x):
    import sys
    sys.path.append('/global/homes/h/huikong/eboss/LSSanalysis')
    from Cosmo import distance
    d = distance(om,1.-om)
    D = d.D(z)
    f = d.omz(z)**.557
    xi = np.loadtxt('./xi0isoChallenge_matterpower6.0.dat').transpose()[1]
    y = np.zeros(len(x))
    for i in range(0,len(x)):
        for j in range(0,len(xi0)):
           if x[i] <= xi0[j]:
              if j == 0:
                  y[i] = xi[j]*x[i]/xi0[j]
                  y[i] = y[i]*x[i]*x[i]
                  break
              else:
                  y[i] = xi[j-1]*(xi0[j]-x[i])/(xi0[j]-xi0[j-1])+xi[j]*(x[i]-xi0[j-1])/(xi0[j]-xi0[j-1])
                  y[i] = y[i]*x[i]*x[i]
                  break
    return y

def see_by_eye():
    b=1
    xi = xigg(b)
    plt.plot(xi0,xi*xi0*xi0)
    a =  np.loadtxt('./NewData_subfiles_chunk21.txt').transpose()
    plt.errorbar(a[0],a[1],yerr = a[2])
    plt.show()
    b = input('b = ?')
    while(b!=0):
       plt.clf()
       xi = xigg(b)
       plt.plot(xi0,xi*xi0*xi0)
       plt.errorbar(a[0],a[1],yerr = a[2])
       plt.show()
       b = input('b = ?')
    plt.clf()
    b = input('save b:\n') 
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages('galaxy_bias_chunk21.pdf') as pdf:
        plt.title('galaxy bias: '+str(b))
        xi = xigg(b)
        plt.plot(xi0,xi*xi0*xi0)
        plt.errorbar(a[0],a[1],yerr = a[2])
        pdf.savefig()
        plt.close()
    return True

def curve_fit():
   from scipy.optimize import curve_fit
   xi = np.loadtxt('./xi0isoChallenge_matterpower6.0.dat').transpose()[1]
   ydata = np.loadtxt(InputFile).transpose()[1]
   sigma_data = np.loadtxt(InputFile).transpose()[2]
   popt, pcov = curve_fit(xigg,xdata,ydata,sigma = sigma_data)
   print('popt = '+str(popt))
   from Cosmo import distance
   d = distance(om,1.-om)
   D = d.D(z)
   print 'D = '+str(D)
   f = d.omz(z)**.557
   p1 = plt.plot(xi0, xi*xi0*xi0,color = 'red')
   xi = xi/(1.+2/3.*.4+0.2*.4**2.)
   xi = popt**2.*(1.+2/3.*f/popt + 0.2*(f/popt)**2)*xi*D
   p2 = plt.plot(xi0, xi*xi0*xi0,color = 'green')
   import matplotlib.patches as mpatches
   red_patch = mpatches.Patch(color='red', label='no bias')
   green_patch = mpatches.Patch(color='green', label='bias:'+str(popt[0]))
   blue_patch  = mpatches.Patch(color='blue', label= 'Chunk22 Data')
   plt.legend(handles=[red_patch,green_patch,blue_patch])
   plt.errorbar(xdata,ydata,sigma_data,color = 'blue')
   from matplotlib.backends.backend_pdf import PdfPages
   with PdfPages('./output/'+OutputFile) as pdf:
       plt.title('galaxy bias: '+str(popt[0])+' +/-'+str(pcov[0][0]))
       plt.xlabel('Mpc')
       pdf.savefig()   
       plt.show()

def chi2(b,d,e,indmin,indmax):
        y = xigg(xdata,b)
        mod = y[indmin:indmax]
        dmod = mod-d
        chi2 = sum((dmod/e)**2.)
        return chi2
        
def findchi2min(rmin=20,rmax=200,bmin=0.9,bmax=1.5,bstep=.001):     
        data = np.loadtxt(InputFile).transpose()
        r = data[0]
        dmod = []
        imin = 0
        imax = 0
        sm = 0
        sx = 0
        for i in range(0,len(r)):
                if r[i] > rmin and r[i] < rmax:
                        if sm == 0:
                                imin = i
                                sm = 1
                if r[i] > rmax and sx == 0:
                        imax = i
                        sx = 1
        if sx == 0:
           imax = len(r)-1
           sx=1
        b = bmin
        d = data[1][imin:imax]
        e = data[2][imin:imax]
        c2min = 1000000
        while b < bmax:
                c2 = chi2(b,d,e,imin,imax)
                if c2 < c2min:
                        c2min = c2
                        bmin = b
                print b,c2           
                b += bstep  
        sup_b = bmin
        inf_b = bmin
        sup_c2 = c2min
        inf_c2 = c2min
        while sup_c2<1+c2min:
            sup_b += bstep
            sup_c2 = chi2(sup_b,d,e,imin,imax)
        print 'sup_b,sup_c2'
        print sup_b,sup_c2
        while inf_c2<1+c2min:
             inf_b -= bstep
             inf_c2 = chi2(inf_b,d,e,imin,imax)
        print 'inf_b,inf_c2'
        print inf_b,inf_c2
        b_err_ave = (sup_b-inf_b)/2
        from Cosmo import distance
        d = distance(om,1.-om)
        D = d.D(z)
        print 'D = '+str(D)
        f = d.omz(z)**.557
        #xi = np.loadtxt('./xi0isoChallenge_matterpower6.0.dat').transpose()[1]
        #p0 = plt.plot(xi0,xi*xi0*xi0)
        p1 = plt.plot(xdata, xigg2(xdata),'r.')
        #xi = xi/(1.+2/3.*.4+0.2*.4**2.)
        #xi = bmin**2.*(1.+2/3.*f/bmin + 0.2*(f/bmin)**2)*xi*D
        p2 = plt.plot(xdata, xigg(xdata,bmin),'g.')
        import matplotlib.patches as mpatches
        red_patch = mpatches.Patch(color='red', label='no bias')
        green_patch = mpatches.Patch(color='green', label='bias:'+str(bmin)+r'$\pm$'+str(b_err_ave))
        blue_patch  = mpatches.Patch(color='blue', label= 'Chunk22 Data')
        plt.legend(handles=[red_patch,green_patch,blue_patch])
        ydata = np.loadtxt(InputFile).transpose()[1]
        sigma_data = np.loadtxt(InputFile).transpose()[2]
        plt.errorbar(xdata,ydata,sigma_data,color = 'blue')
        plt.errorbar
        plt.axis([0,200,-50,100])
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages('./output/'+OutputFile) as pdf:
            plt.title('galaxy bias: '+str(bmin)+r'$\pm$'+str(b_err_ave))
            plt.xlabel('Mpc')
            pdf.savefig()
            plt.show()

def findchi2min_part(Filename,rmin=20,rmax=200,bmin=0.9,bmax=1.5,bstep=.001):
        data = np.loadtxt(Filename).transpose()
        r = data[0]
        dmod = []
        imin = 0
        imax = 0
        sm = 0
        sx = 0
        for i in range(0,len(r)):
                if r[i] > rmin and r[i] < rmax:
                        if sm == 0:
                                imin = i
                                sm = 1
                if r[i] > rmax and sx == 0:
                        imax = i
                        sx = 1
        if sx == 0:
           imax = len(r)-1
           sx=1
        b = bmin
        d = data[1][imin:imax]
        e = data[2][imin:imax]
        c2min = 1000000
        while b < bmax:
                c2 = chi2(b,d,e,imin,imax)
                if c2 < c2min:
                        c2min = c2
                        bmin = b
                #print b,c2           
                b += bstep  
        sup_b = bmin
        inf_b = bmin
        sup_c2 = c2min
        inf_c2 = c2min
        while sup_c2<1+c2min:
            sup_b += bstep
            sup_c2 = chi2(sup_b,d,e,imin,imax)
        print 'sup_b,sup_c2'
        print sup_b,sup_c2
        while inf_c2<1+c2min:
             inf_b -= bstep
             inf_c2 = chi2(inf_b,d,e,imin,imax)
        print 'inf_b,inf_c2'
        print inf_b,inf_c2
        b_err_ave = (sup_b-inf_b)/2
        return bmin,b_err_ave
        #from Cosmo import distance
        #d = distance(om,1.-om)
        #D = d.D(z)
        #print 'D = '+str(D)
        #f = d.omz(z)**.557
        #xi = np.loadtxt('./xi0isoChallenge_matterpower6.0.dat').transpose()[1]
        #p0 = plt.plot(xi0,xi*xi0*xi0)
def upweightdata():
        bmin,b_err_ave = findchi2min_part(InputFile)
        bmin2,b_err_ave2 = findchi2min_part(InputFile2)
        p1 = plt.plot(xdata, xigg(xdata,1),'r.')
        #xi = xi/(1.+2/3.*.4+0.2*.4**2.)
        #xi = bmin**2.*(1.+2/3.*f/bmin + 0.2*(f/bmin)**2)*xi*D
        p2 = plt.plot(xdata, xigg(xdata,bmin),'g.')
        p21= plt.plot(xdata,xigg(xdata,bmin2),'y.')
        import matplotlib.patches as mpatches
        red_patch = mpatches.Patch(color='red', label='no bias')
        green_patch = mpatches.Patch(color='green', label='up weight bias:'+str(bmin)+r'$\pm$'+str(b_err_ave))
        yellow_patch = mpatches.Patch(color='yellow',label='down weight bias:'+str(bmin2)+r'$\pm$'+str(b_err_ave2))
        blue_patch  = mpatches.Patch(color='blue', label= 'up weight ELG data')
        cyan_patch = mpatches.Patch(color='cyan', label= 'down weight ELG random')
        plt.legend(handles=[red_patch,green_patch,yellow_patch,blue_patch,cyan_patch])
        ydata = np.loadtxt(InputFile).transpose()[1]
        sigma_data = np.loadtxt(InputFile).transpose()[2]
        plt.plot(xdata,ydata,color = 'blue')
        plt.fill_between(xdata,ydata-sigma_data,ydata+sigma_data,color='blue',alpha=0.5)
        ydata2 = np.loadtxt(InputFile2).transpose()[1]
        sigma_data2 = np.loadtxt(InputFile2).transpose()[2]
        plt.plot(xdata,ydata,color = 'cyan')
        plt.fill_between(xdata,ydata2-sigma_data2,ydata2+sigma_data2,color='cyan',alpha=0.5)
        plt.axis([0,200,-50,100])
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages('./output/'+OutputFile) as pdf:
            plt.title('Combined Correlation Function and Galaxy Bias for 3 chunks')
            plt.xlabel('Mpc')
            pdf.savefig()
            plt.show()
            
if __name__=='__main__' :
      upweightdata()
