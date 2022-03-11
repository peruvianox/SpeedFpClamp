"""
SpeedFpClamp Data Analysis 
Ricky Pimentel
UNC Chapel Hill Applied Biomechanics Laboratory
2022

Run the script to perform data analysis and generate all article figures
Data avaible at https://drive.google.com/drive/folders/1-u74AgFj0rZDlK9dtHf0nTUbUMXmFZLc?usp=sharing
"""

import pandas as pd
import numpy as np
import os
import fnmatch
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime as dt
import scipy.io as sio
import scipy.stats as stats
import pingouin as pg

# identify folder with data files - RENAME to your folder path!
folder = r'D:\UNC_ABL\FpMetabolics_2020\MetData'
# folder = r'D:\Fp_Metabolics\MetData'
files = os.listdir(folder)
os.chdir(folder)
sty = 'seaborn'
mpl.style.use(sty)

#%% RMR Analysis
# initialize dict that will hold metabolic data
Subjects = {
    's001': {'RMR Avg':'NaN'},
    's002': {'RMR Avg':'NaN'},
    's003': {'RMR Avg':'NaN'},
    's004': {'RMR Avg':'NaN'},
    's005': {'RMR Avg':'NaN'},
    's006': {'RMR Avg':'NaN'},
    's007': {'RMR Avg':'NaN'},
    's008': {'RMR Avg':'NaN'},
    's009': {'RMR Avg':'NaN'},
    's010': {'RMR Avg':'NaN'},
    's011': {'RMR Avg':'NaN'},
    's012': {'RMR Avg':'NaN'},
    's013': {'RMR Avg':'NaN'},
    's014': {'RMR Avg':'NaN'},
    's015': {'RMR Avg':'NaN'},
    's016': {'RMR Avg':'NaN'},
    's017': {'RMR Avg':'NaN'},
    's018': {'RMR Avg':'NaN'},
    's019': {'RMR Avg':'NaN'},
    's020': {'RMR Avg':'NaN'},
    }

SubjNames = list(Subjects.keys())
print(SubjNames)

TrialAvg = {}
v = plt.get_cmap('jet')
cNorm  = mpl.colors.Normalize(vmin=0, vmax=len(Subjects))
scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap='jet')

# get RMR file 
pattern = '*REE*'
matching = fnmatch.filter(files, pattern)
d = 0
for i in matching:
    for s in SubjNames:
        if s in i:
            Subj = s
            break 
        
    F = folder + '\\' + i
    RMR = pd.read_excel(F)
    
    # pull VO2 and time data
    VO2_kg = RMR.loc[2:len(RMR),'VO2/Kg']
    VO2 = RMR.loc[2:len(RMR),'VO2']
    VCO2 = RMR.loc[2:len(RMR),'VCO2']
    t = RMR.loc[2:len(RMR),'t'].values
    W = (VO2 /1000 * 16.5 + VCO2 /1000 * 4.51) * 1000 / 60 
    
    # find rows after 3 min
    T = []
    c = 0
    for i in t:
        c = c + 1
        if i.minute >=3:
             T.append(c)  
             
    # calculate average RMR and make array
    AvgVO2_kg = np.mean(VO2_kg[T])
    AvgVO2 = np.mean(VO2[T])
    AvgVCO2 = np.mean(VCO2[T])
    AvgW = np.mean(W[T])
    Seq = list(range(0, len(T)))
    A = np.ones((len(T)))
    AvgRMR_array = A * AvgW
    colorVal = scalarMap.to_rgba(d)
    plt.plot(W[T], color=colorVal, lw=2)
    plt.plot(T, AvgRMR_array, color=colorVal, lw=2)
    plt.show()
            
    Subjects[Subj]['RMR VO2_kg Avg'] = AvgVO2_kg
    Subjects[Subj]['RMR VO2 Avg'] = AvgVO2
    Subjects[Subj]['RMR VCO2 Avg'] = AvgVCO2
    Subjects[Subj]['RMR W Avg'] = AvgW

    print(s + ' Average RMR = ')
    print(AvgVO2_kg)
    print('mL/kg/min')
    d = d+1
    
plt.title('Resting Metabolic Rate')
plt.ylabel('W')
plt.xlabel('Record #')
plt.savefig('RMR.jpg', dpi=300)

#%% Load and extract Active VO2 data
AvgS = []
AvgF = []
plt.close('all')
# get active trial file 
pattern = '*CPET*'
matching = fnmatch.filter(files, pattern)
fig = plt.figure(figsize=[12,10])
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
MkrSz = 10

d = 0
for i in matching:
    colorVal = scalarMap.to_rgba(d)
    
    for s in SubjNames:
        if s in i:
            Subj = s
            break 
        
    if '$' in i:
        break
        
    print(s)
    Folder = folder + '\\' + i
    Active = pd.read_excel(Folder)

    # pull VO2 and time data
    VO2_kg = Active.loc[2:len(Active),'VO2/Kg'].values.tolist()
    VO2 = Active.loc[2:len(Active),'VO2'].values.tolist()
    VCO2 = Active.loc[2:len(Active),'VO2'].values.tolist()
    t = Active.loc[2:len(Active),'t'].values.tolist()
    Mkr = Active.loc[2:len(Active),'Marker'].values.tolist()
    
    
    # convert start and end times to datetimes
    today = dt.datetime.today()
    ProtocolStartTime = dt.datetime.combine(today, t[0])
    ProtocolEndTime = dt.datetime.combine(today, t[-1])
    
    # Active Trial Analysis
    def TrlMetAnalysis(TrlStr, t, VO2_kg, VO2, VCO2, Mkr):
        today = dt.datetime.today()
        
        # add start and end times
        StartStr = 'start' + TrlStr
        EndStr = 'end' + TrlStr
            
        # pull start and end times of trial
        StartInd = Mkr.index(StartStr)
        EndInd = Mkr.index(EndStr)
        
        # define VO2 during trial
        TrialVO2_kg = VO2_kg[StartInd:EndInd]
        TrialVO2 = VO2[StartInd:EndInd]
        TrialVCO2 = VCO2[StartInd:EndInd]
        TrialW = (np.multiply(np.divide(TrialVO2, 1000), 16.58) + 
                  np.multiply(np.divide(TrialVCO2, 1000), 4.51)) * 1000 / 60 
        
        # create start and end times as datetimes
        StartTime = dt.datetime.combine(today, t[StartInd])
        # EndTime = dt.datetime.combine(today, t[EndInd])
        #  TrialTime = EndTime-StartTime
        
        # convert to seconds
        TrlSec = []
        for i in t[StartInd:EndInd]:
            TrlTime = dt.datetime.combine(today, i)
            ts = TrlTime - StartTime
            TrlSec.append(ts.total_seconds())
        
        # find final 2 min of trial
        Final2Min = [x for x in TrlSec if x >= 180]
        Final2MinInd = TrlSec.index(Final2Min[0])
        
        # average VO2 over final 2 min 
        TrialVO2_kgAvg = np.mean(TrialVO2_kg[Final2MinInd:EndInd])
        TrialVO2Avg = np.mean(TrialVO2[Final2MinInd:EndInd])
        TrialVCO2Avg = np.mean(TrialVCO2[Final2MinInd:EndInd])
        TrialWAvg = np.mean(TrialW[Final2MinInd:EndInd])
        
        VO2Data = {
            'Trial Name' : TrlStr,
            'VO2_kg Avg' : TrialVO2_kgAvg,
            'All VO2_kg Data' : TrialVO2,
            'VO2 Avg' : TrialVO2,
            'All VO2 Data' : TrialVO2Avg,
            'VCO2 Avg' : TrialVCO2Avg,
            'All VCO2 Data' : TrialVCO2,
            'Time Values' : TrlSec,
            'W Avg' : TrialWAvg,
            'All TrialW Data' : TrialW,
            }
        
        return VO2Data
    
    # extract & analyze rest times   
    StartInds = []
    i = 0
    for v in Mkr:
        if 'start' in str(v):
            StartInds.append(i)
        i = i + 1
        
    EndInds = []
    i = 0
    for v in Mkr:
        if 'end' in str(v):
            EndInds.append(i)
        i = i + 1
        
    del(StartInds[0])
    del(EndInds[-1])
    RestTime = []
    for i in range(len(StartInds)):
        starts = dt.datetime.combine(today, t[StartInds[i]])
        ends = dt.datetime.combine(today, t[EndInds[i]])
        NumSec = starts - ends
        RestTime.append(NumSec.seconds)

    Subjects[s]['RestTime'] = RestTime
    
    
    # analyze each active trial individually
    S_M20 = TrlMetAnalysis('S_M20', t, VO2_kg, VO2, VCO2, Mkr)
    S_M10 = TrlMetAnalysis('S_M10', t, VO2_kg, VO2, VCO2, Mkr)
    S_Norm = TrlMetAnalysis('S_Norm', t, VO2_kg, VO2, VCO2, Mkr)
    S_P10 = TrlMetAnalysis('S_P10', t, VO2_kg, VO2, VCO2, Mkr)
    S_P20 = TrlMetAnalysis('S_P20', t, VO2_kg, VO2, VCO2, Mkr)
    
    F_M20 = TrlMetAnalysis('F_M20', t, VO2_kg, VO2, VCO2, Mkr)
    F_M10 = TrlMetAnalysis('F_M10', t, VO2_kg, VO2, VCO2, Mkr)
    F_Norm = TrlMetAnalysis('F_Norm', t, VO2_kg, VO2, VCO2, Mkr)
    F_P10 = TrlMetAnalysis('F_P10', t, VO2_kg, VO2, VCO2, Mkr)
    F_P20 = TrlMetAnalysis('F_P20', t, VO2_kg, VO2, VCO2, Mkr)
        
    # get subject mass & Height
    Subjects[s]['Mass'] = Active.loc[5, 'Unnamed: 1']
    Subjects[s]['Height'] = Active.loc[4, 'Unnamed: 1'] / 100
    
    # create arrays with all 5 trials
    # net VO2/kg 
    Y_S = [S_M20['VO2_kg Avg'], 
           S_M10['VO2_kg Avg'], 
           S_Norm['VO2_kg Avg'],
           S_P10['VO2_kg Avg'], 
           S_P20['VO2_kg Avg']]
    Y_S_net = Y_S - Subjects[Subj]['RMR VO2_kg Avg']
    Y_F = [F_M20['VO2_kg Avg'], 
           F_M10['VO2_kg Avg'], 
           F_Norm['VO2_kg Avg'],
           F_P10['VO2_kg Avg'], 
           F_P20['VO2_kg Avg'] ]
    Y_F_net = Y_F - Subjects[Subj]['RMR VO2_kg Avg']
    
    # net watts
    TrialW_S = [S_M20['W Avg'],
                S_M10['W Avg'], 
                S_Norm['W Avg'],
                S_P10['W Avg'], 
                S_P20['W Avg']]
    TrialW_S_net = (TrialW_S - Subjects[Subj]['RMR W Avg']) / Subjects[s]['Mass']
    TrialW_F = [F_M20['W Avg'], 
                F_M10['W Avg'], 
                F_Norm['W Avg'],
                F_P10['W Avg'], 
                F_P20['W Avg']]
    TrialW_F_net = (TrialW_F - Subjects[Subj]['RMR W Avg']) / Subjects[s]['Mass']
    

    
    # save in dict
    Subjects[s]['Trial_S_VO2'] = Y_S
    Subjects[s]['Trial_S_VO2net'] = Y_S_net 
    Subjects[s]['Trial_F_VO2'] = Y_F
    Subjects[s]['Trial_F_VO2net'] = Y_F_net
    Subjects[s]['TrialW_S_net'] = TrialW_S_net
    Subjects[s]['TrialW_F_net'] = TrialW_F_net

    # plot VO2 data by condition
    X = [1, 2, 3, 4, 5]
    # plt.plot([x - 0.1 for x in X] ,Y_S_net, '.', color=colorVal, lw=5, ms=MkrSz, alpha=0.6,
    #          label = Subj + ' Speed')
    # plt.plot([x + 0.1 for x in X], Y_F_net, '^', color=colorVal, lw=5, ms=MkrSz, alpha=0.6,
    #          label = Subj + ' Force')
    
    # save subject data to aggregate later
    TrialAvg['Subj'] = Subj
    AvgS.append(TrialW_S_net)
    AvgF.append(TrialW_F_net)

    d = d+1
    
# change shape of output into NxCondition numpy array
WAvg_S = np.reshape(AvgS, [len(Subjects), 5])
WAvg_F = np.reshape(AvgF, [len(Subjects), 5])

# calculate averages
TrialAvg['Trial_S_Wnet'] = np.mean(AvgF,axis=0)
TrialAvg['Trial_S_Wnet_sd'] = np.std(AvgF,axis=0)
TrialAvg['Trial_F_Wnet'] = np.mean(AvgF,axis=0)
TrialAvg['Trial_F_Wnet_sd'] = np.std(AvgF,axis=0)

# boxplot
c = 'red'
box = plt.boxplot(WAvg_S, positions=[x - 0.1 for x in X], 
            widths=0.16, patch_artist=True,
            boxprops=dict(facecolor=c, color=c),
            capprops=dict(color=c),
            whiskerprops=dict(color=c),
            flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=c),
            )

c = 'c'
plt.boxplot(WAvg_F, positions=[x + 0.1 for x in X], 
            widths=0.16, patch_artist=True,
            boxprops=dict(color=c, facecolor=c),
            capprops=dict(color=c),
            whiskerprops=dict(color=c),
            flierprops=dict(color=c, markeredgecolor=c),
            medianprops=dict(color=c),
            )

# add labels and such to plot
ax.set_title('Metabolic Cost Across All Trials', fontsize=20)
ax.set_xticks([1, 2, 3, 4, 5])
ax.set_xticklabels(['-20','-10','Norm','+10', '+20'], fontsize=15)
plt.text(1,25, '')
plt.xlabel('Trial', fontsize=15)
plt.ylabel('Net W/kg', fontsize=15)   
plt.show()
# plt.savefig('Normalized Trial VO2.jpg', dpi=300)

#%% Rest Time Analysis
RestTimes = []
for s in Subjects:
    for x in Subjects[s]['RestTime']:
        RestTimes.append(x) 
        
AvgRestTimes = np.mean(RestTimes)
SDRestTimes = np.std(RestTimes)
MinRestTimes = np.min(RestTimes)
MaxRestTimes = np.max(RestTimes)

print('Rest Times')
print('Avg = ' + str(AvgRestTimes))
print('SD = ' + str(SDRestTimes))
print('Min = ' + str(MinRestTimes))
print('Max = ' + str(MaxRestTimes))

#%% Load Matlab Data
plt.close('all')
SubjNamez = []
for v in SubjNames:
    SubjNamez.append(v.replace('s','Subj'))

Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
SubjData = {}

pattern = '*.mat'
matching = fnmatch.filter(files, pattern)
for i in matching:
    for s in SubjNamez:
        if s in i:
            Subj = s
            break 
        
    print('Loading' + s)
    F = folder + '\\' + i
    Dict = {}
    MAT = sio.loadmat(F, mdict = Dict, squeeze_me=1)

    NormSpd = MAT['normSpeed']
    SpdTargets = MAT['speedTargets']
    FpTargets = MAT['FpTargets'][0,:]
    
    def SpdAnalysis(Cond, MAT, Color):
        
        # get variables
        # Data = pd.DataFrame()
        Data = {}
        Spd = MAT['FpTarget'][Cond]['Data'][:]['Speed'].tolist()
        Time = MAT['FpTarget'][Cond]['Data'][:]['Time'].tolist()
        FpTarget = MAT['FpTarget'][Cond]['TargetFp']
        MeanPkFp = MAT['FpTarget'][Cond]['Data'][:]['MeanPeakFp'].tolist()
        MeanPkFb = MAT['NewFpTarget'][Cond]['Data'][:]['MeanPeakFb'].tolist()
        
        # get indicies of final 2 min
        Final2Min = [x for x in Time if x >= 180]
        Final2MinInd = Time.index(Final2Min[Cond])
        EndInd = len(Time)
        
        # calc avg speed over final 2 min
        AvgSpd = np.nanmean(Spd[Final2MinInd:EndInd])
        A = np.ones(len(Time))
        AvgSpd_array = A * AvgSpd
        
        # cacl avg Fp over final 2 min 
        AvgFp = np.nanmean(MeanPkFp[Final2MinInd:EndInd])
        AvgFp_array = A * AvgFp
        AvgFb = np.nanmean(MeanPkFb[Final2MinInd:EndInd])
        AvgFb_array = A * AvgFb
        
        Z = np.zeros(len(Time))
        Z[Final2MinInd:EndInd] = 1
        
        # save data
        Data['Final2min'] = Z
        Data['F_Time'] = Time
        Data['F_Fp'] = MeanPkFp
        Data['F_Fb'] = MeanPkFb
        Data['F_AvgFp'] = AvgFp_array
        Data['F_AvgFb'] = AvgFb_array
        Data['F_Spd'] = Spd
        Data['F_AvgSpd'] = AvgSpd_array
        Data['F_Target'] = FpTarget
        
        
        
        # analyze speed targeting trial  
        Spd = MAT['NewSpeedTarget'][Cond]['Data'][:]['Speed'].tolist()
        Time = MAT['NewSpeedTarget'][Cond]['Data'][:]['Time'].tolist()
        # FpTarget = MAT['NewSpeedTarget'][Cond]['TargetSpeed']
        MeanPkFp = MAT['NewSpeedTarget'][Cond]['Data'][:]['MeanPeakFp'].tolist()
        MeanPkFb = MAT['NewSpeedTarget'][Cond]['Data'][:]['MeanPeakFb'].tolist()
        
        # get indicies of final 2 min
        Final2Min = [x for x in Time if x >= 180]
        Final2MinInd = Time.index(Final2Min[Cond])
        EndInd = len(Time)
        
        # calc avg speed over final 2 min
        AvgSpd = np.nanmean(Spd[Final2MinInd:EndInd])
        A = np.ones(len(Time))
        AvgSpd_array = A * AvgSpd
        
        # cacl avg Fp over final 2 min 
        AvgFp = np.nanmean(MeanPkFp[Final2MinInd:EndInd])
        AvgFp_array = A * AvgFp
        AvgFb = np.nanmean(MeanPkFb[Final2MinInd:EndInd])
        AvgFb_array = A * AvgFb

        Z = np.zeros(len(Time))
        Z[Final2MinInd:EndInd] = 1
        
        # save data
        Data['S_Final2min'] = Z
        Data['S_Time'] = Time
        Data['S_Fp'] = MeanPkFp
        Data['S_Fb'] = MeanPkFb
        Data['S_AvgFp'] = AvgFp_array
        Data['S_AvgFb'] = AvgFb_array
        Data['S_Spd'] = Spd
        Data['S_AvgSpd'] = MAT['SpeedTarget'][Cond]['Speed']
        
        return Data#, SpdData
     
        
    DataM20 = SpdAnalysis(0, MAT, 'blue')
    DataM10 = SpdAnalysis(1, MAT, 'cornflowerblue')
    DataNorm = SpdAnalysis(2, MAT, 'black')
    DataP10 = SpdAnalysis(3, MAT, 'orange')
    DataP20 = SpdAnalysis(4, MAT, 'orangered')
    
    
    SubjData.update({s+'M20': DataM20, 
                s+'M10': DataM10, 
                s+'Norm': DataNorm, 
                s+'P10': DataP10, 
                s+'P20': DataP20})

del DataM20, DataM10, DataNorm, DataP10, DataP20, Dict, MAT, F
del Active, FpTargets, i, EndInds, ends, Mkr
del t, T, today, v, W

#%% Sampling Frequency Analysis
SampFreq = []
for s in SubjData:
    if 'Spd' not in s:
        t = SubjData[s]['F_Time']
        SampFreq.append(np.mean(np.diff(t))) 
    
SamplingMean = np.mean(SampFreq)
SamplingStd = np.std(SampFreq)

print('Avg Sampling Intermission = ' + str(SamplingMean))
print('Std Sampling Intermission = ' + str(SamplingStd))

#%% Plot Fps and Fbs during speed clamp
sty = 'default'
mpl.style.use(sty)
import seaborn as sns
import matplotlib.colors as mcolors
plt.close('all')
fnt = 12
plt.rcParams.update({'font.size': fnt})
fig = plt.figure(figsize=[12,20])
ax0 = plt.subplot(621)
ax1 = plt.subplot(623)
ax2 = plt.subplot(625)
ax3 = plt.subplot(622)
ax4 = plt.subplot(624)
ax5 = plt.subplot(626)
ax6 = plt.subplot(212)

ax0.set_position([0.1, 0.80, 0.36, 0.125])
ax1.set_position([0.1, 0.65, 0.36, 0.125])
ax2.set_position([0.1, 0.45, 0.36, 0.16])
ax3.set_position([0.54, 0.80, 0.36, 0.125])
ax4.set_position([0.54, 0.65, 0.36, 0.125])
ax5.set_position([0.54, 0.45, 0.36, 0.16])
ax6.set_position([0.1, 0.06, 0.8, 0.3])


from matplotlib.patches import Rectangle
Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
Lvl = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['black'],
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
A1 = 0.3
A2 = 1
A3 = 0.05
LW1 = 2
LW2 = 2
counter = 0
SubjFps = np.zeros([len(SubjNamez), 5])
SubjFbs = np.zeros([len(SubjNamez), 5])
SubjFpf = np.zeros([len(SubjNamez), 5])
SubjFbf = np.zeros([len(SubjNamez), 5])
    
for s in SubjNamez:
    for t in [0, 1, 2, 3, 4]:
        # calculate targeting accuracy
        S = s.replace('Subj', 's')
        Weight = Subjects[S]['Mass'] * 9.81
        
        # speed Clamp
        Time = SubjData[s+Levels[t]]['S_Time']
        NormTarget = SubjData[s+'Norm']['F_Target']
        TrlFp = SubjData[s+Levels[t]]['S_Fp']
        TrlFb = SubjData[s+Levels[t]]['S_Fb']
        
        Ind = 200
        Fp = TrlFp[Ind:]
        Fb = TrlFb[Ind:]
        
        ax0.plot(Time[Ind:], 100*np.divide(Fp, Weight),
                  c=Colors[t], alpha=A1)
        ax1.plot(Time[Ind:], 100*np.divide(Fb, Weight), 
          c=Colors[t], alpha=A1)

        shade = mcolors.to_rgb(Colors[t]) + (A3,)
        ax2.plot(100*np.divide(Fp, Weight), 100*np.divide(Fb, Weight), 
                 'o', ms=6, mfc=shade, mec=shade, mew=0)
        
        # get average Fp over final 2 min
        Final2Min = [x for x in Time if x >= 180]
        a = Time.index(Final2Min[0])
        b = len(Time)
        AvgFp = np.ones(len(TrlFp[a:b])) * np.mean(TrlFp[a:b])
        AvgFb = np.ones(len(TrlFb[a:b])) * np.mean(TrlFb[a:b])
        
        SubjFps[counter, t] = 100*AvgFp[0] / Weight
        SubjFbs[counter, t] = 100*AvgFb[0] / Weight
        
        
        # Fp clamp
        Time = SubjData[s+Levels[t]]['F_Time']
        TrlFp = SubjData[s+Levels[t]]['F_Fp']
        TrlFb = SubjData[s+Levels[t]]['F_Fb']
        Ind = 200
        Fp = TrlFp[Ind:]
        Fb = TrlFb[Ind:]
        
        ax3.plot(Time[Ind:], 100*np.divide(Fp, Weight),
                  c=Colors[t], alpha=A1)
        ax4.plot(Time[Ind:], 100*np.divide(Fb, Weight), 
          c=Colors[t], alpha=A1)

    
        shade = mcolors.to_rgb(Colors[t]) + (A3,)
        ax5.plot(100*np.divide(Fp, Weight), 100*np.divide(Fb, Weight), 
                 's', ms=6, mfc=shade, mec=shade, mew=0)
        
        # get average Fp over final 2 min
        Final2Min = [x for x in Time if x >= 180]
        a = Time.index(Final2Min[0])
        # a = Time[Time == Final2Min[0]].index[0]
        b = len(Time)
        AvgFp = np.ones(len(TrlFp[a:b])) * np.mean(TrlFp[a:b])
        AvgFb = np.ones(len(TrlFb[a:b])) * np.mean(TrlFb[a:b])
        
        SubjFpf[counter, t] = 100*AvgFp[0] / Weight
        SubjFbf[counter, t] = 100*AvgFb[0] / Weight
        
        
    counter = counter + 1
    
# for i in range(5):
#     ax0.add_patch(Rectangle(((180,np.mean(SubjFps[:,i])-0.025)), 120, 0.05,
#               edgecolor = Colors[i],
#               facecolor = Colors[i],
#               alpha = A2,
#               fill=True, lw=0))
#     val = int(np.mean(SubjFps[:,i]))
#     font = mpl.font_manager.FontProperties()
#     font.set_weight('bold')
#     font.set_size(fnt)
#     ax0.text(240, np.mean(SubjFps[:,i]), 
#               Lvl[i] + ' Avg: ' + str(val) + '%', 
#               va='center', ha='center', c='w', 
#               fontproperties=font)
    
#     ax1.add_patch(Rectangle(((180,np.mean(SubjFbs[:,i])-0.025)), 120, 0.05,
#               edgecolor = Colors[i],
#               facecolor = Colors[i],
#               alpha = A2,
#               fill=True, lw=0))
#     val = int(np.mean(SubjFbs[:,i]))
#     font = mpl.font_manager.FontProperties()
#     font.set_weight('bold')
#     font.set_size(fnt)
#     ax1.text(240, np.mean(SubjFbs[:,i]), 
#               Lvl[i] + ' Avg: ' + str(val) + '%', 
#               va='center', ha='center', c='w', 
#               fontproperties=font)

g = 0.5
z = np.linspace(0,50)
ax2.plot(z, z, '-', lw=2, c=[g, g, g])
ax5.plot(z, z, '-', lw=2, c=[g, g, g])

# plt.axvline(x = 180, color='k', lw=2)
# plt.text(240, 0.62, 'Final 2 Minutes', c='k', 
#          fontsize=fnt, ha='center')
ax0.set_xlim([0, 300])
ax0.set_ylim([0, 50])
ax0.set_xticks([0, 100, 200, 300])
ax1.set_xlim([0, 300])
ax1.set_ylim([0, 50])
ax1.set_xticks([0, 100, 200, 300])
ax2.set_xlim([0, 50])
ax2.set_ylim([0, 50])
ax3.set_xlim([0, 300])
ax3.set_ylim([0, 50])
ax3.set_xticks([0, 100, 200, 300])
ax4.set_xlim([0, 300])
ax4.set_ylim([0, 50])
ax4.set_xticks([0, 100, 200, 300])
ax5.set_xlim([0, 50])
ax5.set_ylim([0, 50])


# ax0.set_xlabel('Time(s)')
ax0.set_ylabel('Fp (%BW)')
ax0.set_title('Speed Clamp: All Steps')
# ax1.set_title('Speed Clamp Fb')
ax1.set_ylabel('Fb (%BW)')
ax1.set_xlabel('Time(s)')
ax2.set_xlabel('Fp (%BW)')
ax2.set_ylabel('Fb (%BW)')
# ax2.set_title('Speed Clamp Fp : Fb')

# ax3.set_xlabel('Time(s)')
# ax3.set_ylabel('% Body Weight')
ax3.set_title('Fp Clamp: All Steps')
# ax4.set_title('Fp Clamp Fb')
ax4.set_xlabel('Time(s)')

ax5.set_xlabel('Fp (%BW)')
ax5.set_ylabel('Fb (%BW)')
# ax5.set_title('Fp Clamp Fp : Fb')
# plt.title('Fp during Speed Clamp')
# plt.text(-25, 1.42, 'A', fontsize=30)


A1 = 0.5
A2 = 0.5
MS = 8
for i in range(5): 
    ax6.plot(SubjFps[:,i], SubjFbs[:,i], 'o',
             c = Colors[i], alpha=A1, ms=MS)
    ax6.plot(SubjFpf[:,i], SubjFbf[:,i], 's',
             c = Colors[i], alpha=A1, ms=MS)

vals = np.linspace(10, 40)
ax6.plot(vals, vals, '-', c=[g, g, g], lw=2)
ax6.text(35, 22, 'Unity Line', color = [g, g, g], 
         va='center', ha='right', fontsize=fnt+3)
ax6.plot([31, 32], [22, 22], '-', c = [g, g, g])

FpS = SubjFps.reshape((1,100)).squeeze()
FbS = SubjFbs.reshape((1,100)).squeeze()
FpF = SubjFpf.reshape((1,100)).squeeze()
FbF = SubjFbf.reshape((1,100)).squeeze()
ax6.set_xlim((12, 36))
ax6.set_ylim((12, 36))
ax6.set_xlabel('Fp (%BW)')
ax6.set_ylabel('Fb (%BW)')
ax6.set_title('Average Fb & Fp Over Final 2 Minutes')


y = sns.regplot(FpS, FbS, color='k', scatter=False, 
            line_kws={'linestyle':'-'}, ax=ax6, fit_reg=True)
sns.regplot(FpF, FbF, color='k', scatter=False,  
            line_kws={'linestyle':'--'}, ax=ax6)

fnt = 15
from scipy import stats
from scipy.stats import t

# calculate linear regression
res = stats.linregress(FpS, FbS)
tinv = lambda p, df: abs(t.ppf(p/2, df))
ts = tinv(0.05, len(FpS)-2)
d = ts*res.stderr
ax6.text(13, 34, 'Speed Clamp', fontsize=fnt+3, va='center')
ax6.plot(18, 34, 'ok', ms=15)
ax6.plot([18.5, 20], [34, 34], '-k', ms=16)
s = round(res.slope, 3)
l = round(s - d,3)
u = round(s + d,3)
S = 'Slope: ' + str(s) + ' 95% CI: (' + str(l) + ', ' + str(u) + ')'
ax6.text(13, 32, S, fontsize=fnt, va='center')

# calculate linear regression
res = stats.linregress(FpF, FbF)
tinv = lambda p, df: abs(t.ppf(p/2, df))
ts = tinv(0.05, len(FpS)-2)
d = ts*res.stderr
ax6.text(26, 16, 'Fp Clamp', fontsize=fnt+3, va='center')
ax6.plot(30, 16, 'sk', ms=16)
ax6.plot([30.5, 32], [16, 16], '--k', ms=15)
s = round(res.slope, 3)
l = round(s - d,3)
u = round(s + d,3)
S = 'Slope: ' + str(s) + ' 95% CI: (' + str(l) + ', ' + str(u) + ')'
ax6.text(26, 14, S, fontsize=fnt, va='center')


plt.savefig('Fp-Fb.tiff', dpi=300)
# plt.savefig('Fp-Fb.pdf', dpi=300)


#%% Plot Fps during Speed Clamp
sty = 'default'
mpl.style.use(sty)
plt.close('all')
fnt = 16
plt.rcParams.update({'font.size': fnt})
fig = plt.figure(figsize=[20,10])
# ax = plt.subplot(131)
ax = plt.axes((0.06, 0.075, 0.27, 0.75))
from matplotlib.patches import Rectangle
Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
Lvl = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['black'],
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
A1 = 0.15
A2 = 0.7
LW1 = 2
LW2 = 2
counter = 0
SubjFps = np.zeros([len(SubjNamez), 5])
    
for s in SubjNamez:
    for t in [0, 1, 2, 3, 4]:
    # calculate targeting accuracy
        # Time = SubjData[s+'Spd'+Levels[t]]['S_Time']
        # NormTarget = SubjData[s+'Norm']['F_Target'][2]
        # TrlFp = SubjData[s+'Spd'+Levels[t]]['S_Fp']
        Time = SubjData[s+Levels[t]]['S_Time']
        NormTarget = SubjData[s+'Norm']['F_Target']
        TrlFp = SubjData[s+Levels[t]]['S_Fp']
        Fp = list(filter(None, TrlFp))
        Ind = TrlFp.index(Fp[0])
        
        plt.plot(Time[Ind:], np.divide(Fp,NormTarget), 
                  c=Colors[t], alpha=A1)
        
         # get average Fp over final 2 min
        Final2Min = [x for x in Time if x >= 180]
        a = Time.index(Final2Min[0])
        b = len(Time)
        AvgFp = np.ones(len(TrlFp[a:b])) * np.mean(TrlFp[a:b])
        
        SubjFps[counter, t] = AvgFp[0] / NormTarget
    counter = counter + 1
    
for i in range(5):
    ax.add_patch(Rectangle(((180,np.mean(SubjFps[:,i])-0.025)), 120, 0.05,
              edgecolor = 'k',
              facecolor = Colors[i],
              alpha = A2,
              zorder=3,
              fill=True, lw=2))
    val = int(np.mean(SubjFps[:,i])*100)
    font = mpl.font_manager.FontProperties()
    font.set_weight('bold')
    font.set_size(fnt)
    plt.text(240, np.mean(SubjFps[:,i]), 
              Lvl[i] + ' Avg: ' + str(val) + '%', 
              va='center', ha='center', c='w', 
              fontproperties=font)
    
plt.axvline(x = 180, color='k', lw=2)
plt.text(240, 0.62, 'Final 2 Minutes', c='k', 
         fontsize=fnt, ha='center')
plt.xlim([0, 300])
plt.ylim([0.6, 1.4])
ax.set_xticks([0, 100, 200, 300])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4])
ax.set_yticklabels(['60%', '70%', '80%', '90%', '100%',
                  '110%', '120%', '130%', '140%'])

plt.xlabel('Time(s)')
plt.ylabel('Relative Fp')
plt.title('Fp during Speed Clamp')
plt.text(-25, 1.42, 'A', fontsize=30)

#%% Plot Fps during Fp Clamp
# ax = plt.subplot(132)
ax = plt.axes((0.39, 0.075, 0.27, 0.75))
from matplotlib.patches import Rectangle
Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
Lvl = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['black'],
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
A1 = 0.15
A2 = 1
LW1 = 2
LW2 = 2
counter = 0
SubjFps = np.zeros([len(SubjNamez), 5])
def MovingAvg(Vals, Window):
    Filt = Vals
    Win = int((Window-1) / 2)
    for x in np.arange(Win,len(Vals)-Win):
        Filt[x] = np.nanmean([Vals[x-Win:x+Win]])
        
    return Filt
    
for s in SubjNamez:
    for t in [0, 1, 2, 3, 4]:
    # calculate targeting accuracy
        # Time = SubjData[s+Levels[t]]['F_Time'].to_list()
        # NormTarget = SubjData[s+'Norm']['F_Target'][2]
        # TrlFp = SubjData[s+Levels[t]]['F_Fp'].values.tolist() 
        Time = SubjData[s+Levels[t]]['F_Time']
        NormTarget = SubjData[s+'Norm']['F_Target']
        TrlFp = SubjData[s+Levels[t]]['F_Fp']
        Fp = list(filter(None, TrlFp))
        Ind = TrlFp.index(Fp[0])
       
        plt.plot(Time[Ind:],  np.divide(Fp,NormTarget), 
                  c=Colors[t], alpha=A1)
        
         # get average Fp over final 2 min
        Final2Min = [x for x in Time if x >= 180]
        a = Time.index(Final2Min[0])
        b = len(Time)
        AvgFp = np.ones(len(TrlFp[a:b])) * np.mean(TrlFp[a:b])
        
        SubjFps[counter, t] = AvgFp[0] / NormTarget
    counter = counter + 1
    
for i in range(5):
    ax.add_patch(Rectangle(((180,np.mean(SubjFps[:,i])-0.025)), 120, 0.05,
              edgecolor = 'k',
              facecolor = Colors[i],
              alpha = A2,
              zorder=3,
              fill=True, lw=2))
    val = int(np.mean(SubjFps[:,i])*100)
    font = mpl.font_manager.FontProperties()
    font.set_weight('bold')
    font.set_size(fnt)
    plt.text(240, np.mean(SubjFps[:,i]), 
              Lvl[i] + ' Avg: ' + str(val) + '%', 
              va='center', ha='center', c='w', 
              fontproperties=font)
    
plt.axvline(x = 180, color='k', lw=2)
plt.text(240, 0.62, 'Final 2 Minutes', c='k', 
         fontsize=fnt, ha='center')
plt.xlim([0, 300])
plt.ylim([0.6, 1.4])
ax.set_xticks([0, 100, 200, 300])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4])
ax.set_yticklabels(['60%', '70%', '80%', '90%', '100%',
                  '110%', '120%', '130%', '140%'])

plt.xlabel('Time(s)')
plt.ylabel('Relative Fp')
plt.title('Fp during Fp Clamp')
plt.text(-25, 1.42, 'B', fontsize=30)

#%% Plot General Speeds
# from matplotlib.patches import Rectangle
Levels = ['M20', 'M10', 'Norm', 'P10', 'P20']
Lvl = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['black'],
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]

ax = plt.axes((0.72, 0.075, 0.27, 0.75))
A1 = 0.3
A2 = 1
LW1 = 2
LW2 = 2
SubjSpds = np.zeros([len(SubjNamez), 5])
counter = 0
for s in SubjNamez:
    
    NormSpeed = SubjData[s+'Norm']['S_AvgSpd']
    # time = SubjData[s+'Norm']['F_Time'].to_list()
    NormSpd = np.divide(SubjData[s+'Norm']['F_Spd'], NormSpeed)
    
    for t in [0, 1, 2, 3, 4]:
    
        # plot individual speed lines for final 4 min
        time = SubjData[s+Levels[t]]['F_Time']
        TrlSpd = np.divide(SubjData[s+Levels[t]]['F_Spd'], NormSpeed)

        plt.plot(time, TrlSpd, lw=LW1,
                 c=Colors[t], alpha=A1)
        
        # plot final 2 min average speeds
        Final2Min = [x for x in time if x >= 180]
        a = time.index(Final2Min[0])
        b = len(time)
        AvgSpd = np.ones(len(TrlSpd[a:b])) * np.nanmean(TrlSpd[a:b])
        SubjSpds[counter, t] = AvgSpd[0]
    counter = counter + 1

for i in range(5):
    ax.add_patch(Rectangle(((180,np.mean(SubjSpds[:,i])-0.025)), 120, 0.05,
             edgecolor = 'k',
             facecolor = Colors[i],
             alpha = A2,
             zorder=3,
             fill=True, lw=2))
    val = int(np.nanmean(SubjSpds[:,i])*100)
    font = mpl.font_manager.FontProperties()
    font.set_weight('bold')
    font.set_size(fnt)
    plt.text(240, np.mean(SubjSpds[:,i]), 
             Lvl[i] + ' Avg: ' + str(val) + '%', 
             va='center', ha='center', c='w', 
             fontproperties=font)
    
plt.axvline(x = 180, color='k', lw=2)
plt.text(240, 0.62, 'Final 2 Minutes', c='k', 
         fontsize=fnt, ha='center')
plt.xlim([0, 300])
plt.ylim([0.6, 1.4])
ax.set_xticks([0, 100, 200, 300])
ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4])
ax.set_yticklabels(['60%', '70%', '80%', '90%', '100%',
                  '110%', '120%', '130%', '140%'])

plt.xlabel('Time(s)')
plt.ylabel('Relative Walking Speed')
plt.text(-25, 1.42, 'C', fontsize=30)
plt.title('Walking Speed during Fp Clamp')

plt.savefig('BiofeedbackPerformance.tiff', dpi=300)
# plt.savefig('BiofeedbackPerformance.pdf', dpi=300)
    

#%% plot Speed and Fp from each trial 
plt.close('all')
fig = plt.figure(figsize=[12,12])
Conditions = ['M20', 'M10', 'P10', 'P20']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
TrialInd = [0, 1, 3, 4]
AllSpd_F = []
AllFp_F = []
AllSpd_S = []
AllFp_S = []
AllFb_S = []
AllFb_F = []
N = len(Subjects)
Ones = np.ones(N)
Mass = [0]*N
MassTxt = [0]*N
Alls = [0]*N
C = 0
for s in Subjects:
    subjSpd_S = [0, 0, 0, 0, 0]
    subjFp_S = [0, 0, 0, 0, 0]
    subjSpd_F = [0, 0, 0, 0, 0]
    subjFp_F = [0, 0, 0, 0, 0]
    subjFb_F = [0, 0, 0, 0, 0]
    subjFb_S = [0, 0, 0, 0, 0]
    
    Mass[C] = Subjects[s]['Mass']
    MassTxt[C] = s

    Subj = s.replace('s', 'Subj')
    Trial = 'Norm'
    Key1 = Subj + Trial
    # Key2 = Subj + Trial
    
    # generate norms
    NormSpd_F = SubjData[Key1]['F_AvgSpd'][0]
    NormFp_F = SubjData[Key1]['F_AvgFp'][0]
    NormSpd_S = SubjData[Key1]['S_AvgSpd']
    NormFp_S = SubjData[Key1]['S_AvgFp'][0]
    NormFb_S = SubjData[Key1]['S_AvgFb'][0]
    NormFb_F = SubjData[Key1]['F_AvgFb'][0]
    
    subjSpd_S[2] = NormSpd_S
    subjFp_S[2] = NormFp_S
    subjSpd_F[2] = NormSpd_F
    subjFp_F[2] = NormFp_F
    subjFb_F[2] = NormFb_F
    subjFb_S[2] = NormFb_S
    
    ax1 = fig.add_subplot(221)
    plt.scatter(1, Subjects[s]['TrialW_S_net'][2],
                    c='k', marker='.')
    ax2 = fig.add_subplot(222)
    plt.scatter(1, Subjects[s]['TrialW_S_net'][2],
                    c='k', marker='.')
    
    ax3 = fig.add_subplot(223)
    plt.scatter(1, Subjects[s]['TrialW_F_net'][2],
                    c='k', marker='.')
    ax4 = fig.add_subplot(224)
    plt.scatter(1, Subjects[s]['TrialW_F_net'][2],
                    c='k', marker='.')
    
    # loop through non-norm conditions
    for cond in [0, 1, 2, 3]:
    # for cond in [0, 1, 3, 4]:
        Trial = Conditions[cond]
        Key1 = Subj + Trial
        # Key2 = Subj + Trial
        
        # calculate & normalize variables
        Spd_S = SubjData[Key1]['S_AvgSpd']
        Fp_S = SubjData[Key1]['S_AvgFp'][0]
        Spd_F = SubjData[Key1]['F_AvgSpd'][0]
        Fp_F = SubjData[Key1]['F_AvgFp'][0]
        Fb_F = SubjData[Key1]['F_AvgFb'][0]
        Fb_S = SubjData[Key1]['S_AvgFb'][0]

        subjSpd_S[TrialInd[cond]] = Spd_S
        subjFp_S[TrialInd[cond]] = Fp_S
        subjSpd_F[TrialInd[cond]] = Spd_F
        subjFp_F[TrialInd[cond]] = Fp_F
        subjFb_F[TrialInd[cond]] = Fb_F
        subjFb_S[TrialInd[cond]] = Fb_S
        
          # plot values being sure to normalize to Norm trial
        ax1 = fig.add_subplot(221)
        plt.scatter(Spd_S/NormSpd_S,
                    Subjects[s]['TrialW_S_net'][TrialInd[cond]], 
                    c=Colors[cond], marker='.', alpha=0.5)
        ax2 = fig.add_subplot(222)
        plt.scatter(Fp_S/NormFp_S,
                    Subjects[s]['TrialW_S_net'][TrialInd[cond]], 
                    c=Colors[cond], marker='.', alpha=0.5)
        

        ax3 = fig.add_subplot(223)
        plt.scatter(Spd_F/NormSpd_F, 
                        Subjects[s]['TrialW_F_net'][TrialInd[cond]], 
                        c=Colors[cond], marker='.', alpha=0.5)
        ax4 = fig.add_subplot(224)
        plt.scatter(Fp_F/NormFp_F, 
                        Subjects[s]['TrialW_F_net'][TrialInd[cond]], 
                        c=Colors[cond], marker='.', alpha=0.5)
        
        
    AllSpd_S.append(subjSpd_S)
    AllFp_S.append(subjFp_S)
    AllSpd_F.append(subjSpd_F)
    AllFp_F.append(subjFp_F)
    AllFb_F.append(subjFb_F)
    AllFb_S.append(subjFb_S)
    Alls[C] = Subj
    C = C+1
    
    
AllSpd_S = np.reshape(AllSpd_S, [len(Subjects), 5])
AllFp_S = np.reshape(AllFp_S, [len(Subjects), 5])
AllSpd_F = np.reshape(AllSpd_F, [len(Subjects), 5])
AllFp_F = np.reshape(AllFp_F, [len(Subjects), 5])
AllFb_F = np.reshape(AllFb_F, [len(Subjects), 5])
AllFb_S = np.reshape(AllFb_S, [len(Subjects), 5])

CoT_Fp_S = WAvg_S / AllSpd_S
CoT_Fp_F = WAvg_F / AllSpd_F
Fp_S = np.array(AllFp_S)
Fp_F = np.array(AllFp_F)
Fb_S = np.array(AllFb_S)
Fb_F = np.array(AllFb_F)

for i in range(len(Subjects)):
    Fp_S[i,:] = AllFp_S[i,:] / Mass[i]
    Fp_F[i,:] = AllFp_F[i,:] / Mass[i]
    Fb_S[i,:] = AllFb_S[i,:] / Mass[i]
    Fb_F[i,:] = AllFb_F[i,:] / Mass[i]

Fp_S = np.reshape(Fp_S, [len(Subjects), 5])
Fp_F = np.reshape(Fp_F, [len(Subjects), 5])
Fb_S = np.reshape(Fb_S, [len(Subjects), 5])
Fb_F = np.reshape(Fb_F, [len(Subjects), 5])

    
ax1.set_xlabel('Normalized Speed', fontsize=10)
ax1.set_ylabel('Net W/kg', fontsize=10)
ax1.set_title('Speed during Fixed Speed', fontsize=12)

ax2.set_xlabel('Normalized Fp', fontsize=10)
# ax2.set_ylabel(ax, 'Net W/kg', fontsize=10)
ax2.set_title('Fp during Fixed Speed', fontsize=12)

ax3.set_xlabel('Normalized Speed', fontsize=10)
ax3.set_ylabel('Net W/kg', fontsize=10)
ax3.set_title('Speed during Fp Targeting', fontsize=12)

ax4.set_xlabel('Normalized Fp', fontsize=10)
# ax4.set_ylabel(ax, 'Net W/kg', fontsize=10)
ax4.set_title('Fp during Fp Targeting', fontsize=12)

# plt.savefig('SpeedsFps.jpg', dpi=300)


#%% Speed Between Fixed and Targeting
plt.close('all')
Conditions = ['-20', '-10', 'Norm', '+10', '+20']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          'k',
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
N = len(Subjects)
Ones = np.ones(N)
Mass = [0]*N
Ofst = 0.1
BarOfst = 0.2
Trans = 0.4
Trans2 = 1
MkrSz = 16
Fnt = 12
TFnt = 16

fig = plt.figure(figsize=[12,12])
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
AllFp_S_kg = np.zeros_like(AllSpd_S)
AllFp_F_kg = np.zeros_like(AllSpd_S)

AllFb_S_kg = np.zeros_like(AllSpd_S)
AllFb_F_kg = np.zeros_like(AllSpd_S)

for x in range(5):
    for i in range(N):
        ax1.plot([X[x]-Ofst, X[x]+Ofst], 
                  [AllSpd_S[i, x], AllSpd_F[i,x]],
                  '-', c=Colors[x], alpha=Trans)
        
        Mass[i] = Subjects[SubjNames[i]]['Mass']
        AllFp_S_kg[i, x] = 9.81 * AllFp_S[i, x] / Mass[i]
        AllFp_F_kg[i, x] = 9.81 * AllFp_F[i, x] / Mass[i]
        
        AllFb_S_kg[i, x] = 9.81 * AllFb_S[i, x] / Mass[i]
        AllFb_F_kg[i, x] = 9.81 * AllFb_F[i, x] / Mass[i]
        
        ax2.plot([X[x]-Ofst, X[x]+Ofst], 
                  [AllFp_S[i, x] / Mass[i], AllFp_F[i,x] / Mass[i]],
                  '-', c=Colors[x], alpha=Trans)
        
        ax3.plot(AllSpd_S[i, x], AllSpd_F[i,x], 
                  '.', c=Colors[x], alpha=Trans2)
        
        ax4.plot(AllFp_S[i, x] / Mass[i], AllFp_F[i,x] / Mass[i], 
                  '.', c=Colors[x], alpha=Trans2)
        
    # plot group averages
    ax1.errorbar(X[x]-BarOfst, np.mean(AllSpd_S[:, x], axis=0), 
                  yerr=np.std(AllSpd_S[:, x], axis=0),      
                  marker='.', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    ax1.errorbar(X[x]+BarOfst, np.mean(AllSpd_F[:, x], axis=0), 
                  yerr=np.std(AllSpd_F[:, x], axis=0),      
                  marker='^', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    ax2.errorbar(X[x]-BarOfst, np.mean(AllFp_S_kg[:, x], axis=0), 
                  yerr=np.std(AllFp_S_kg[:, x], axis=0),      
                  marker='.', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    ax2.errorbar(X[x]+BarOfst, np.mean(AllFp_F_kg[:, x], axis=0), 
                  yerr=np.std(AllFp_F_kg[:, x], axis=0),      
                  marker='^', c=Colors[x], ecolor=Colors[x], markersize=MkrSz)
    
ax1.set_xlabel('Condition', fontsize=Fnt)
ax1.set_ylabel('m/s', fontsize=Fnt)
ax1.set_title('Speed Across Conditions', fontsize=TFnt)
ax1.set_xticks(X)
ax1.set_xticklabels(Conditions)

ax2.set_xlabel('Condition', fontsize=Fnt)
ax2.set_ylabel('N / kg', fontsize=Fnt)
ax2.set_title('Fp Across Conditions', fontsize=TFnt)
ax2.set_xticks(X)
ax2.set_xticklabels(Conditions)

ax3.set_xlabel('Fixed Speed (m/s)', fontsize=Fnt)
ax3.set_ylabel('Fp Targeting (m/s)', fontsize=Fnt)
ax3.set_title('Speed by Condition', fontsize=TFnt)

ax4.set_xlabel('Fixed Speed (N/kg)', fontsize=Fnt)
ax4.set_ylabel('Fp Targeting (N/kg)', fontsize=Fnt)
ax4.set_title('Fp by Condition', fontsize=TFnt)

# plt.savefig('SpeedsFpComp.jpg', dpi=300)


#%% Abstract Plot (Speed, Fp, and CoT)
plt.close('all')
Conditions = ['-20%', '-10%', 'Norm', '+10%', '+20%']
Colors = [mpl._color_data.CSS4_COLORS['blue'], 
          mpl._color_data.CSS4_COLORS['cornflowerblue'], 
          'k',
          mpl._color_data.CSS4_COLORS['orange'], 
          mpl._color_data.CSS4_COLORS['orangered']]
N = len(Subjects)
Ones = np.ones(N)
# Mass = [0]*N
Ofst = 0.12
Trans = 0.25
Full = 1
MkrSz = 10
MkrSz2 = 14
fnt = 15
X = [1, 2, 3, 4, 5]

fig = plt.figure(figsize=[12,12])
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

for x in range(5):
        
    # 1st plot - speed
    ax1.plot([X[x]-Ofst, X[x]+Ofst], 
             [AllSpd_S[:, x], AllSpd_F[:,x]],
             '-', c=Colors[x], alpha=Trans)
    ax1.errorbar(X[x]-BarOfst, np.mean(AllSpd_S[:, x], axis=0), 
                 yerr=np.std(AllSpd_S[:, x], axis=0),      
                 marker='o', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    ax1.errorbar(X[x]+BarOfst, np.mean(AllSpd_F[:, x], axis=0), 
                 yerr=np.std(AllSpd_F[:, x], axis=0),      
                 marker='s', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    
    # 2nd plot - Fp
    ax2.plot([X[x]-Ofst, X[x]+Ofst], 
             [AllFp_S_kg[:, x], AllFp_F_kg[:,x]],
             '-', c=Colors[x], alpha=Trans)
    ax2.errorbar(X[x]-BarOfst, np.mean(AllFp_S_kg[:, x], axis=0), 
                 yerr=np.std(AllFp_S_kg[:, x], axis=0),      
                 marker='o', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    ax2.errorbar(X[x]+BarOfst, np.mean(AllFp_F_kg[:, x], axis=0), 
                 yerr=np.std(AllFp_F_kg[:, x], axis=0),      
                 marker='s', c=Colors[x], ecolor=Colors[x], markersize=MkrSz2)
    
    # 3rd plot - Met Cost
    ax3.plot([X[x]-Ofst, X[x]+Ofst], 
             [WAvg_S[:,x], WAvg_F[:,x]],
             '-', c=Colors[x], alpha=Trans)
    ax3.errorbar(X[x]-BarOfst, np.mean(WAvg_S[:,x], axis=0), 
         yerr=np.std(WAvg_S[:,x], axis=0), 
         marker='o', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    ax3.errorbar(X[x]+BarOfst, np.mean(WAvg_F[:,x], axis=0), 
         yerr=np.std(WAvg_F[:,x], axis=0), 
         marker='s', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    
    # 4th plot - CoT
    ax4.plot([X[x]-Ofst, X[x]+Ofst], 
             [CoT_Fp_S[:,x], CoT_Fp_F[:,x]],
             '-', c=Colors[x], alpha=Trans)
    ax4.errorbar(X[x]-BarOfst, np.mean(CoT_Fp_S[:,x], axis=0), 
         yerr=np.std(CoT_Fp_S[:,x], axis=0), 
         marker='o', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    ax4.errorbar(X[x]+BarOfst, np.mean(CoT_Fp_F[:,x], axis=0), 
         yerr=np.std(CoT_Fp_F[:,x], axis=0), 
         marker='s', ecolor=Colors[x], color=Colors[x], markersize=MkrSz2)
    
    
# create legend
ax1.plot(1, 1.95, 'o', color='k', markersize=MkrSz2)
ax1.plot(1, 1.85, 's', color='k', markersize=MkrSz2)
ax1.text(1.12, 1.95, ' = Speed Clamp', color='k', fontsize=18, va='center')
ax1.text(1.12, 1.85, ' = Fp Clamp', color='k', fontsize=18, va='center')
ax1.text(0, 2.08, 'A', fontsize=fnt*2)
ax2.text(0, 30.6, 'B', fontsize=fnt*2)
ax3.text(0, 10.5, 'C', fontsize=fnt*2)
ax4.text(0, 6.2, 'D', fontsize=fnt*2)

# edit axes
# ax1.set_xlabel('Condition', fontsize=15)
ax1.set_ylabel('Walking Speed (m/s)', fontsize=15)
# ax1.set_title('A', fontsize=20, horizontalalignment='left')
ax1.set_xticks(X)
ax1.set_xticklabels(Conditions, fontsize=15)
ax1.tick_params(axis='y', labelsize=15) 
# plt.title(label='A', fontsize=20, Loc='left')
# ax1.set_ylim(1, 2)

# ax2.set_xlabel('Condition', fontsize=15)
ax2.set_ylabel('Fp (%BW)', fontsize=15)
# ax2.set_title('B', fontsize=20, horizontalalignment='left')
ax2.set_xticks(X)
ax2.set_xticklabels(Conditions, fontsize=15)
ax2.tick_params(axis='y', labelsize=15) 
ax2.set_ylim(14, 30)


# ax3.plot(1, 9.2, 'o', color='k', markersize=MkrSz2)
# ax3.plot(1, 8.2, 's', color='k', markersize=MkrSz2)
# ax3.text(1.12, 9.2, ' = Speed Clamp', color='k', fontsize=18, va='center')
# ax3.text(1.12, 8.2, ' = Fp Clamp', color='k', fontsize=18, va='center')

ax3.set_xticks(X)
ax3.set_xticklabels(Conditions, fontsize=15)
ax3.set_ylabel('Net Metabolic Power (W/kg)', fontsize=15)
# ax3.set_title('C', fontsize=20, Loc='left')
ax3.tick_params(axis='y', labelsize=15) 
ax3.tick_params(axis='x', labelsize=15)
# ax3.set_xlim(1.3, 3)
ax3.set_ylim(0, 10)

ax4.set_xticks(X)
ax4.set_xticklabels(Conditions, fontsize=15)
ax4.set_ylabel('CoT (J/kg/m)', fontsize=15)
# ax4.set_title('D', fontsize=20, Loc='left')
ax4.tick_params(axis='y', labelsize=15) 
ax4.set_ylim(1, 6)
ax4.tick_params(axis='x', labelsize=15)


#%% Run stats and add to plot
S = range(1,len(Subjects)+1)
Ones = np.ones(5)
fnt = 12

#perform the repeated measures ANOVA
RMA = pd.DataFrame({'subjects': np.tile(np.repeat(S, len(X)), 2),
                   'condition': np.tile(X, len(Subjects)*2),
                   'clamp': np.repeat(np.hstack((Ones, Ones*2)), len(Subjects)),
                   'speed': np.reshape([AllSpd_S, AllSpd_F],
                                       [len(Subjects)*2*5, 1][0]),
                   'Fp': np.reshape([AllFp_S_kg, AllFp_F_kg],
                                    [len(Subjects)*2*5, 1][0]), 
                   'MetCost': np.reshape([WAvg_S, WAvg_F],
                                         [len(Subjects)*2*5, 1][0]), 
                   'CoT': np.reshape([CoT_Fp_S, CoT_Fp_F],
                                     [len(Subjects)*2*5, 1][0])}
                    )

AnovaNames = ['speed', 'Fp', 'MetCost', 'CoT']
aov = {}
for A in AnovaNames:
    aov[A] = pg.rm_anova(data=RMA, dv=A, within=['condition', 'clamp'], subject='subjects', detailed=True)
    print('\n\n' + A + '\n')
    print('P values: ')
    print(aov[A]['p-unc'])
    print('Partial Eta Sq: ')
    print(aov[A]['np2'])
    


# place ANOVA values in fig
ax1.text(3.8, 1.15,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
         va='top', fontsize = fnt, ha='right')
ax1.text(4.0, 1.19,'     p'+'\n'+'<0.001'+'\n'+'  0.004'+'\n'+'<0.001', 
         va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['speed']['np2'].to_list(), 3)
ax1.text(4.75, 1.225,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')


ax2.text(3.8, 16.5,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
         va='top', fontsize = fnt, ha='right')
ax2.text(4, 17,'     p'+'\n'+'<0.001'+'\n'+'  0.001'+'\n'+'<0.001', 
         va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['Fp']['np2'].to_list(), 3)
ax2.text(4.75, 17,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')

ax3.text(3.8, 1.675,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
          va='top', fontsize = fnt, ha='right')
ax3.text(4, 2,'     p'+'\n'+'<0.001'+'\n'+'  0.002'+'\n'+'  0.126', 
          va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['MetCost']['np2'].to_list(), 3)
ax3.text(4.75, 2.325,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')

ax4.text(3.8, 1.8,'condition'+'\n'+'clamp'+'\n'+'condition x clamp', 
         va='top', fontsize = fnt, ha='right')
ax4.text(4, 2,'     p'+'\n'+'<0.001'+'\n'+'  0.010'+'\n'+'  0.313', 
         va='top', fontsize = fnt, ha='left')
np2 = np.round(aov['CoT']['np2'].to_list(), 3)
ax4.text(4.75, 2.175,'    $\eta^2_p$\n'+str(np2[0])+'\n'+str(np2[1])+'\n'+str(np2[2]), 
         va='top', fontsize = fnt, ha='left')

#%% Post hoc T-tests
Ast = 26
Hash = 18

# speed sub-analysis
# test across conditions
T_SCondSpeed = np.ones(5)
T_FCondSpeed = np.ones(5)
ES_SCondSpeed = np.ones(5)
ES_FCondSpeed = np.ones(5)
G = np.array([np.ones(20), 2*np.ones(20)]).reshape([40, 1])
for x in [0, 1, 3, 4]:
    d = np.reshape([AllSpd_S[:,x], AllSpd_S[:,2]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_SCondSpeed[x] = float(Stats['p-tukey'])
    ES_SCondSpeed[x] = float(Stats['cohen'])
    if T_SCondSpeed[x] < 0.05 :
        ax1.text(x+1-BarOfst, np.mean(AllSpd_S[:,x])+0.18, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')
        
    d = np.reshape([AllSpd_F[:,x], AllSpd_F[:,2]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_FCondSpeed[x] = float(Stats['p-tukey'])
    ES_FCondSpeed[x] = float(Stats['cohen'])
    if T_FCondSpeed[x] < 0.05 :
        ax1.text(x+1+BarOfst, np.mean(AllSpd_F[:,x])+0.18, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')
        
    
print('\nSpeed Speed Conditions Post Hoc')
print('p-values')
print(T_SCondSpeed)
print('effect sizes')
print(np.round(ES_SCondSpeed, decimals=5))
print('Speed Fp Conditions Post Hoc')
print('p-values')
print(np.round(T_FCondSpeed, decimals=5))
print('effect sizes')
print(np.round(ES_FCondSpeed, decimals=5))

# between clamps
T_Speed = np.ones(5)
ES_Speed = np.ones(5)
for x in range(5):
    d = np.reshape([AllSpd_S[:,x], AllSpd_F[:,x]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_Speed[x] = float(Stats['p-tukey'])
    ES_Speed[x] = float(Stats['cohen'])
    if T_Speed[x] < 0.05 :
        y = np.mean([np.mean(AllSpd_S[:,x], axis=0), np.mean(AllSpd_F[:,x], axis=0)])
        ax1.text(x+1, y+0.2, '#', 
                 c = Colors[x], fontsize=Hash, ha='center')

print('Speed Between Clamp Post Hoc')
print('p-values')
print(np.round(T_Speed, decimals=5))
print('effect sizes')
print(np.round(ES_Speed, decimals=5))
print(' ')
print(' ')
        
# Fp sub-analysis
# between conditions
T_SCondFp = np.ones(5)
T_FCondFp = np.ones(5)
ES_SCondFp = np.ones(5)
ES_FCondFp = np.ones(5)
for x in [0, 1, 3, 4]:
    d = np.reshape([AllFp_S_kg[:,x], AllFp_S_kg[:,2]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_SCondFp[x] = float(Stats['p-tukey'])
    ES_SCondFp[x] = float(Stats['cohen'])
    if T_SCondFp[x] < 0.05 :
        ax2.text(x+1-BarOfst, np.mean(AllFp_S_kg[:,x])+2.5, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')

    d = np.reshape([AllFp_F_kg[:,x], AllFp_F_kg[:,2]], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_FCondFp[x] = float(Stats['p-tukey'])
    ES_FCondFp[x] = float(Stats['cohen'])
    if T_FCondFp[x] < 0.05 :
        ax2.text(x+1+BarOfst, np.mean(AllFp_F_kg[:,x], axis=0)+2.5, '*', 
                 c = Colors[x], fontsize=Ast, ha='center')
        
print('Fp Speed Conditions Post Hoc')
print('p-values')
print(T_SCondFp)
print('effect sizes')
print(np.round(ES_SCondFp, decimals=5))
print('Fp Fp Conditions Post Hoc')
print('p-values')
print(T_FCondFp)
print('effect sizes')
print(np.round(ES_FCondFp, decimals=5))

# between clamps
T_Fp = np.ones(5)
ES_Fp = np.ones(5)
for x in range(5):
    d = np.reshape([AllFp_S_kg[:,x], AllFp_F_kg[:,x]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_Fp[x] = float(Stats['p-tukey'])
    ES_Fp[x] = float(Stats['cohen'])
    if T_Fp[x] < 0.05 :
        y = np.mean([np.mean(AllFp_S_kg[:,x], axis=0), 
                     np.mean(AllFp_F_kg[:,x], axis=0)])
        ax2.text(x, y+3, '#', c = Colors[x], fontsize=Hash, ha='center')
print('Fp Between Clamp Post Hoc')
print('p-values')
print(np.round(T_Speed, decimals=5))
print('effect sizes')
print(np.round(ES_Speed, decimals=5))
print(' ')
print(' ')

# MetCost sub-analysis
# between conditions
T_SCondMetCost = np.ones(5)
T_FCondMetCost = np.ones(5)
ES_SCondMetCost = np.ones(5)
ES_FCondMetCost = np.ones(5)
for x in range(5): 
    d = np.reshape([WAvg_S[:,x], WAvg_S[:,2]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_SCondMetCost[x] = float(Stats['p-tukey'])
    ES_SCondMetCost[x] = float(Stats['cohen'])
    if T_SCondMetCost[x] < 0.05 :
        ax3.text(x+1-BarOfst, np.mean(WAvg_S[:,x], axis=0)+1.5,
                 '*', c = Colors[x], fontsize=Ast, ha='center')
    
    d = np.reshape([WAvg_F[:,x], WAvg_F[:,2]], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_FCondMetCost[x] = float(Stats['p-tukey'])
    ES_FCondMetCost[x] = float(Stats['cohen'])
    if T_FCondMetCost[x] < 0.05 :
        ax3.text(x+1+BarOfst, np.mean(WAvg_F[:,x], axis=0)+1.5, 
                 '*', c = Colors[x], fontsize=Ast, ha='center')
        
print('MetCost Speed Conditions Post Hoc')
print('p-values')
print(T_SCondMetCost)
print('effect sizes')
print(np.round(ES_SCondMetCost, decimals=5))
print('MetCost Conditions Post Hoc')
print('p-values')
print(T_FCondMetCost)
print('effect sizes')
print(np.round(ES_FCondMetCost, decimals=5))

# post hoc difference in net metabolic cost for lowest condition intensity
# np.mean([WAvg_F[:,0]]) - np.mean([WAvg_S[:,0]])

# between clamps
T_MetCost = np.ones(5)
ES_MetCost = np.ones(5)
for x in range(5):
    d = np.reshape([WAvg_S[:,x], WAvg_F[:,x]], [40, 1])    
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_MetCost[x] = float(Stats['p-tukey'])
    ES_MetCost[x] = float(Stats['cohen'])
    if T_MetCost[x] < 0.05 :
        y = np.mean([np.mean(WAvg_S[:,x], axis=0), 
                     np.mean(WAvg_F[:,x], axis=0)])
        ax3.text(x+1, y+1.7, '#', 
                 c = Colors[x], fontsize=Hash, ha='center')

print('MetCost Between Clamp Post Hoc')
print('p-values')
print(np.round(T_MetCost, decimals=5))
print('effect sizes')
print(np.round(ES_MetCost, decimals=5))
print(' ')
print(' ')

# CoT sub-analysis
# between conditions
T_SCondCoT = np.ones(5)
T_FCondCoT = np.ones(5)
ES_SCondCoT = np.ones(5)
ES_FCondCoT = np.ones(5)
for x in range(5):
    d = np.reshape([CoT_Fp_S[:,x], CoT_Fp_S[:,2]], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_SCondCoT[x] = float(Stats['p-tukey'])
    ES_SCondCoT[x] = float(Stats['cohen'])
    if T_SCondCoT[x] < 0.05 :
        ax4.text(x+1-BarOfst, np.mean(CoT_Fp_S[:,x], axis=0)+1,
                 '*', c = Colors[x], fontsize=Ast, ha='center')
    
    d = np.reshape([CoT_Fp_F[:,x], CoT_Fp_F[:,2]], [40, 1]) 
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_FCondCoT[x] = float(Stats['p-tukey'])
    ES_FCondCoT[x] = float(Stats['cohen'])
    if T_FCondCoT[x] < 0.05 :
        ax4.text(x+1+BarOfst, np.mean(CoT_Fp_F[:,x], axis=0)+1, 
                 '*', c = Colors[x], fontsize=Ast, ha='center')
        
print('CoT Speed Conditions Post Hoc')
print('p-values')
print(T_SCondCoT)
print('effect sizes')
print(np.round(ES_SCondCoT, decimals=5))
print('CoT Fp Conditions Post Hoc')
print('p-values')
print(T_FCondCoT)
print('effect sizes')
print(np.round(ES_FCondCoT, decimals=5))

# between clamps
T_CoT = np.ones(5)
ES_CoT = np.ones(5)
for x in range(5):
    d = np.reshape([CoT_Fp_S[:,x], CoT_Fp_F[:,x]], [40, 1])
    D = pd.DataFrame(np.hstack([d,G]), columns=['X','G'])
    Stats = pg.pairwise_tukey(D, dv='X', between='G', effsize='cohen')
    T_CoT[x] = float(Stats['p-tukey'])
    ES_CoT[x] = float(Stats['cohen'])
    if T_CoT[x] < 0.05 :
        y = np.mean([np.mean(CoT_Fp_S[:,x], axis=0), 
                     np.mean(CoT_Fp_F[:,x], axis=0)])
        ax4.text(x+1, y+1.1, '#', 
                 c = Colors[x], fontsize=Hash, ha='center')

print('CoT Between Clamp Post Hoc')
print('p-values')
print(np.round(T_CoT, decimals=5))
print('effect sizes')
print(np.round(ES_CoT, decimals=5))
print(' ')
print(' ')   

plt.savefig('Clamps.tiff', dpi=300)
# plt.savefig('Clamps.pdf', dpi=300)

#%% quantify differences between clamp types
Fd = np.zeros(5)
Sd = np.zeros(5)
NMPd = np.zeros(5)
CoTd = np.zeros(5)
for x in range(5):
    Sd[x] = abs(np.mean(AllSpd_S[:,x]) - np.mean(AllSpd_F[:,x])) / (0.5 * (np.mean(AllSpd_S[:,x]) + np.mean(AllSpd_F[:,x])))
    Fd[x] = abs(np.mean(AllFp_S_kg[:,x]) - np.mean(AllFp_F_kg[:,x])) / (0.5 * (np.mean(AllFp_S_kg[:,x]) + np.mean(AllFp_F_kg[:,x])))             
    NMPd[x] = abs(np.mean(WAvg_S[:,x]) - np.mean(WAvg_F[:,x])) / (0.5 * (np.mean(WAvg_S[:,x]) + np.mean(WAvg_F[:,x])))
    CoTd[x] = abs(np.mean(CoT_Fp_S[:,x]) - np.mean(CoT_Fp_F[:,x])) / (0.5 * (np.mean(CoT_Fp_S[:,x]) + np.mean(CoT_Fp_F[:,x])))                 

print('Mean Difference for walking speed = ' + str(np.mean(Sd)))
print('Mean Difference for Fp = ' + str(np.mean(Fd)))
print('Mean Difference for net matabolic power = ' + str(np.mean(NMPd)))
print('Mean Difference for CoT = ' + str(np.mean(CoTd)))


#%% Walking Speed Predictors Correlation Plot
plt.close('all')
fig = plt.figure(figsize=[12,12])
sz = 50
sz2 = 100
A = 0.4
fnt = 15
txt = 13

# define correlation plotting function
def CorrPlot(X1, Y1, X2, Y2, ax, Xlim, Ylim, Invert):
    # 1 = speed clamp, 2 = Fp clamp
    for i in range(5): # scatter plot all measures
        ax.scatter(X1[:,i], Y1[:,i], 
                    c=Colors[i], marker='o', s=sz)
        ax.scatter(X2[:,i], Y1[:,i], 
                    c=Colors[i], marker='s', s=sz)
        
    # speed clamp     
    # calculate trendlines
    z_S = np.polyfit(np.hstack(X1), np.hstack(Y1), 1)
    p_S = np.poly1d(z_S)
    x = np.linspace(np.min(X1), np.max(X1), 25)
    ax.scatter(x,p_S(x),c='k',marker='o', s=sz2, alpha = A)
    ax.plot(x,p_S(x),'-k')
    # the line equation and R
    [c_S, P_S] = stats.pearsonr(np.hstack(X1), np.hstack(Y1))
    LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x + ' + str(round(z_S[1],2))
    if P_S < 0.001:
        t = ' Speed Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3)) + ' \n p < 0.001 '
    else: 
        t = ' Speed Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3))
    if Invert == 'Yes':
        ax.text(Xlim[1], Ylim[1], t, fontsize=txt, ha='right', va='top')
    else:
        ax.text(Xlim[1], Ylim[0], t, fontsize=txt, ha='right', va='bottom')

    
    # Fp clamp
    z_F = np.polyfit(np.hstack(X2), np.hstack(Y2), 1)
    p_F = np.poly1d(z_F)
    x = np.linspace(np.min(X2), np.max(X2), 25)
    ax.scatter(x,p_F(x),c='k',marker='s', s=sz2, alpha = A)
    ax.plot(x,p_F(x),'-k')
    # the line equation and R
    [c_F, P_F] = stats.pearsonr(np.hstack(X2), np.hstack(Y2))
    LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x + ' + str(round(z_F[1],2))
    if P_F < 0.001:
        t = ' Fp Clamp \n ' + LineEq_F + ' \n R$^2$ = ' + str(round(c_F*c_F,3)) + ' \n p < 0.001 '
    else: 
        t = ' Fp Clamp \n ' + LineEq_F + ' \n R$^2$ = ' + str(round(c_F*c_F,3))
    if Invert == 'Yes':
        ax.text(Xlim[0], Ylim[0], t, fontsize=txt, ha='left', va='bottom')
    else:
        ax.text(Xlim[0], Ylim[1], t, fontsize=txt, ha='left', va='top')
    
    ax.set_xlim(Xlim)
    ax.set_ylim(Ylim)


# speed by Fp
plt1 = plt.subplot(221)
CorrPlot(AllSpd_S, AllFp_S_kg, AllSpd_F, AllFp_F_kg,
         plt1, (0.95, 2.05), (12, 34), 'N')

plt1.text(0.85, 34.5, 'A', fontsize=fnt*2)
plt1.set_xlabel('Speed (m/s)', fontsize=fnt)
plt1.set_ylabel('Fp (%BW)', fontsize=fnt)
plt1.set_yticks([15, 20, 25, 30])



# speed by Fb
plt2 = plt.subplot(222)
CorrPlot(AllSpd_S, AllFb_S_kg, AllSpd_F, AllFb_F_kg,
         plt2, (0.95, 2.05), (12, 34), 'N')

plt2.text(0.85, 34.5, 'B', fontsize=fnt*2)
plt2.set_xlabel('Speed (m/s)', fontsize=fnt)
plt2.set_xticks([1, 1.5, 2])
plt2.set_ylabel('Fb (%BW)', fontsize=fnt)
plt2.set_yticks([15, 20, 25, 30])


# speed by Stride Length
SL = pd.read_csv(r'C:\Users\richa\Documents\Bioengineering\Abstracts and Manuscripts\Pimentel SpeedFpClamp\JoB Revision\StrideLength.csv')

# get subject heights
Hts = []
for s in Subjects:    
    Hts.append(Subjects[s]['Height'])

# normalize stride length to % height
SLnorm = np.empty(np.shape(SL))
for i in range(10):
    SLnorm[:,i] = np.divide(SL.iloc[:,i], Hts)
            
plt3 = plt.subplot(223)
CorrPlot(AllSpd_S, SLnorm[:, 5:10], AllSpd_F, SLnorm[:, 0:5],
         plt3, (0.95, 2.05), (0.6, 1), 'N')

plt3.text(0.85, 1.02, 'C', fontsize=fnt*2)
plt3.set_xlabel('Speed (m/s)', fontsize=fnt)
plt3.set_xticks([1, 1.5, 2])
plt3.set_ylabel('Stride Length (% height)', fontsize=fnt)
plt3.set_yticks([0.6, 0.7, 0.8, 0.9, 1])



# speed by Stride Duration
SD = pd.read_csv(r'C:\Users\richa\Documents\Bioengineering\Abstracts and Manuscripts\Pimentel SpeedFpClamp\JoB Revision\StrideDur.csv')
plt4 = plt.subplot(224)
CorrPlot(AllSpd_S, np.array(SD.iloc[:, 5:10]), AllSpd_F, np.array(SD.iloc[:, 0:5]),
         plt4, (0.95, 2.05), (0.75, 1.35), 'Yes')


plt4.text(0.85, 1.375, 'D', fontsize=fnt*2)
plt4.set_xlabel('Speed (m/s)', fontsize=fnt)
plt4.set_xticks([1, 1.5, 2])
plt4.set_ylabel('Stride Duration (s)', fontsize=fnt)
plt4.set_yticks([0.8, 0.9, 1, 1.1, 1.2, 1.3])


plt.savefig('CorrSpeed.tiff', dpi=300)
# plt.savefig('CorrSpeed.pdf', dpi=300)


#%% Correlations to walking economy

# define correlation plotting function
def CorrPlot2(X1, Y1, X2, Y2, ax, Xlim, Ylim, Invert):
    # 1 = speed clamp, 2 = Fp clamp
    for i in range(5): # scatter plot all measures
        ax.scatter(X1[:,i], Y1[:,i], 
                    c=Colors[i], marker='o', s=sz)
        ax.scatter(X2[:,i], Y1[:,i], 
                    c=Colors[i], marker='s', s=sz)
        
    # speed clamp     
    # calculate trendlines
    z_S = np.polyfit(np.hstack(X1), np.hstack(Y1), 2)
    p_S = np.poly1d(z_S)
    x = np.linspace(np.min(X1), np.max(X1), 25)
    ax.scatter(x,p_S(x),c='k',marker='o', s=sz2, alpha = A)
    ax.plot(x,p_S(x),'-k')
    # the line equation and R
    [c_S, P_S] = stats.pearsonr(np.hstack(X1), np.hstack(Y1))
    LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x$^2$ + ' + str(round(z_S[1],2)) + 'x + ' + str(round(z_S[2],2))
    if P_S < 0.001:
        t = ' Speed Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3)) + ' \n p < 0.001 '
    else: 
        t = ' Speed Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3))
    if Invert == 'Yes':
        ax.text(Xlim[1], Ylim[1], t, fontsize=txt, ha='right', va='top')
    else:
        ax.text(Xlim[1], Ylim[0], t, fontsize=txt, ha='right', va='bottom')

    
    # Fp clamp
    z_F = np.polyfit(np.hstack(X2), np.hstack(Y2), 2)
    p_F = np.poly1d(z_F)
    x = np.linspace(np.min(X2), np.max(X2), 25)
    ax.scatter(x,p_F(x),c='k',marker='s', s=sz2, alpha = A)
    ax.plot(x,p_F(x),'-k')
    # the line equation and R
    [c_F, P_F] = stats.pearsonr(np.hstack(X2), np.hstack(Y2))
    LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x$^2$ + ' + str(round(z_F[1],2)) + 'x + ' + str(round(z_F[2],2))
    if P_F < 0.001:
        t = ' Fp Clamp \n ' + LineEq_F + ' \n R$^2$ = ' + str(round(c_F*c_F,3)) + ' \n p < 0.001 '
    else: 
        t = ' Fp Clamp \n ' + LineEq_F + ' \n R$^2$ = ' + str(round(c_F*c_F,3))
    if Invert == 'Yes':
        ax.text(Xlim[0], Ylim[0], t, fontsize=txt, ha='left', va='bottom')
    else:
        ax.text(Xlim[0], Ylim[1], t, fontsize=txt, ha='left', va='top')
    
    ax.set_xlim(Xlim)
    ax.set_ylim(Ylim)



plt.close('all')
fig = plt.figure(figsize=[12,12])


plt1 = plt.subplot(221)
CorrPlot(AllSpd_S, WAvg_S, AllSpd_F, WAvg_F,
         plt1, (0.95, 2.05), (1, 9), 'N')

plt1.text(0.88, 9.3, 'A', fontsize=fnt*2)
plt1.set_xlabel('Speed (m/s)', fontsize=fnt)
plt1.set_xticks([1, 1.5, 2])
plt1.set_ylabel('Net Metabolic Power (W/kg)', fontsize=fnt)
plt1.set_yticks([2, 4, 6, 8])




# speed by CoT
plt2 = plt.subplot(222)
CorrPlot2(AllSpd_S, CoT_Fp_S, AllSpd_F, CoT_Fp_F,
          plt2, (0.95, 2.05), (1, 6), 'N')

plt2.text(0.85, 6.2, 'B', fontsize=fnt*2)
plt2.set_xlabel('Speed (m/s)', fontsize=fnt)
plt2.set_xticks([1, 1.5, 2])
plt2.set_ylabel('Cost of Transport (J/kg/m)', fontsize=fnt)
plt2.set_yticks([2, 3, 4, 5, 6])
plt2.tick_params(axis='y', labelsize=fnt) 


# Fp by net metabolic cost
plt3 = plt.subplot(223)
CorrPlot(AllFp_S_kg, WAvg_S, AllFp_F_kg, WAvg_S,
         plt3, (14, 29), (1, 9.5), 'N')

plt3.text(12.5, 9.7, 'C', fontsize=fnt*2)
plt3.set_xlabel('Fp (%BW)', fontsize=fnt)
plt3.set_xticks([15, 20, 25])
plt3.set_ylabel('Net Metabolic Power (W/kg)', fontsize=fnt)
plt3.set_yticks([2, 4, 6, 8])



# Fp by CoT
plt4 = plt.subplot(224)
CorrPlot2(AllFp_S_kg, CoT_Fp_S, AllFp_F_kg, CoT_Fp_F,
          plt4, (14, 29), (1, 6), 'N')


plt4.text(12.5, 6.2, 'D', fontsize=fnt*2)
plt4.set_xlabel('Fp (%BW)', fontsize=fnt)
plt4.set_xlim([14, 29])
plt4.set_xticks([15, 20, 25])
plt4.set_ylabel('Cost of Transport (J/kg/m)', fontsize=fnt)
plt4.set_yticks([1, 2, 3, 4, 5, 6])
plt4.tick_params(axis='y', labelsize=fnt) 


plt.savefig('CorrEconomy.tiff', dpi=300)
# plt.savefig('CorrEconomy.pdf', dpi=300)

#%% Plot Correlates to stride length and stride duration
fnt = 12
plt.rcParams.update({'font.size': fnt})
plt.close('all')
fig = plt.figure(figsize=[20, 10])
sz = 50
sz2 = 100
A = 0.4
fnt = 14
txt = 11

def CorrPlot(X1, Y1, X2, Y2, ax, Xlim, Ylim, Invert):
    # 1 = speed clamp, 2 = Fp clamp
    for i in range(5): # scatter plot all measures
        ax.scatter(X1[:,i], Y1[:,i], 
                    c=Colors[i], marker='o', s=sz)
        ax.scatter(X2[:,i], Y1[:,i], 
                    c=Colors[i], marker='s', s=sz)
        
    # speed clamp     
    # calculate trendlines
    z_S = np.polyfit(np.hstack(X1), np.hstack(Y1), 1)
    p_S = np.poly1d(z_S)
    x = np.linspace(np.min(X1), np.max(X1), 25)
    ax.scatter(x,p_S(x),c='k',marker='o', s=sz2, alpha = A)
    ax.plot(x,p_S(x),'-k')
    # the line equation and R
    [c_S, P_S] = stats.pearsonr(np.hstack(X1), np.hstack(Y1))
    LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x + ' + str(round(z_S[1],2))
    if P_S < 0.001:
        t = ' Speed Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3)) + ' \n p < 0.001 '
    else: 
        t = ' Speed Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3))
    if Invert == 'Yes':
        ax.text(Xlim[1], Ylim[1], t, fontsize=txt, ha='right', va='top')
    else:
        ax.text(Xlim[1], Ylim[0], t, fontsize=txt, ha='right', va='bottom')

    
    # Fp clamp
    z_F = np.polyfit(np.hstack(X2), np.hstack(Y2), 1)
    p_F = np.poly1d(z_F)
    x = np.linspace(np.min(X2), np.max(X2), 25)
    ax.scatter(x,p_F(x),c='k',marker='s', s=sz2, alpha = A)
    ax.plot(x,p_F(x),'-k')
    # the line equation and R
    [c_F, P_F] = stats.pearsonr(np.hstack(X2), np.hstack(Y2))
    LineEq_F = 'y = ' + str(round(z_F[0],2)) + 'x + ' + str(round(z_F[1],2))
    if P_F < 0.001:
        t = ' Fp Clamp \n ' + LineEq_F + ' \n R$^2$ = ' + str(round(c_F*c_F,3)) + ' \n p < 0.001 '
    else: 
        t = ' Fp Clamp \n ' + LineEq_F + ' \n R$^2$ = ' + str(round(c_F*c_F,3))
    if Invert == 'Yes':
        ax.text(Xlim[0], Ylim[0], t, fontsize=txt, ha='left', va='bottom')
    else:
        ax.text(Xlim[0], Ylim[1], t, fontsize=txt, ha='left', va='top')
    
    ax.set_xlim(Xlim)
    ax.set_ylim(Ylim)
        
        
# Stride Length Correlations
# SL by Fp 
plt1 = plt.subplot(241)
CorrPlot(AllFp_S_kg, SLnorm[:,5:10], AllFp_F_kg, SLnorm[:,0:5], 
         plt1, (14, 30), (0.6, 1), 'N')

plt1.set_ylabel('Stride Length (%height)', fontsize=fnt)
plt1.set_yticks([0.6, 0.7, 0.8, 0.9, 1])

# SL by Fb
plt2 = plt.subplot(242)
CorrPlot(AllFb_S_kg, SLnorm[:,5:10], AllFb_F_kg, SLnorm[:,0:5], 
         plt2, (12, 35), (0.6, 1), 'N')
plt2.set_yticks([0.6, 0.7, 0.8, 0.9, 1])

# SL by NMP
plt3 = plt.subplot(243)
CorrPlot(WAvg_S, SLnorm[:,5:10], WAvg_F, SLnorm[:,0:5],
         plt3, (1, 10), (0.6, 1), 'N')
plt3.set_yticks([0.6, 0.7, 0.8, 0.9, 1])

# SL by CoT
plt4 = plt.subplot(244)
CorrPlot(CoT_Fp_S, SLnorm[:,5:10], CoT_Fp_F, SLnorm[:,0:5],
         plt4, (1, 6), (0.6, 1), 'N')
plt4.set_yticks([0.6, 0.7, 0.8, 0.9, 1])



# Stride Duration
# SD by Fp          
plt5 = plt.subplot(245)   
CorrPlot(AllFp_S_kg, np.array(SD.iloc[:,5:10]), AllFp_F_kg, np.array(SD.iloc[:,0:5]), 
         plt5, (14, 30), (0.8, 1.4), 'Yes')
plt5.set_xlabel('Fp (%BW)', fontsize=fnt)
plt5.set_ylabel('Stride Duration (s)', fontsize=fnt)


# SD by Fb
plt6 = plt.subplot(246)        
CorrPlot(AllFb_S_kg, np.array(SD.iloc[:,5:10]), AllFb_F_kg, np.array(SD.iloc[:,0:5]), 
         plt6, (12, 35), (0.8, 1.4), 'Yes')

plt6.set_xlabel('Fb (%BW)', fontsize=fnt)
# plt6.set_xticks([1, 1.5, 2])

# SD by NMP
plt7 = plt.subplot(247)   
CorrPlot(WAvg_S, np.array(SD.iloc[:,5:10]), WAvg_F, np.array(SD.iloc[:,0:5]), 
         plt7, (1, 10), (0.8, 1.4), 'Yes')
plt7.set_xlabel('Net Metabolic Power (W/kg)', fontsize=fnt)

# SD by CoT
plt8 = plt.subplot(248)        
CorrPlot(CoT_Fp_S, np.array(SD.iloc[:,5:10]), CoT_Fp_F, np.array(SD.iloc[:,0:5]), 
         plt8, (1, 6), (0.8, 1.4), 'Yes')
plt8.set_xlabel('CoT (J/kgm)', fontsize=fnt)


fig.savefig('CorrSDSL.tiff', dpi=300)
# plt.savefig('CorrSDSL.pdf', dpi=300)

#%% Variation in speed correlation with metabolic cost
plt.close('all')
plt.figure(figsize = (16, 8))
ax = plt.subplot(121)
ax2 = plt.subplot(122)
SpeedVar = np.zeros((20,5))
# get speed variability per trial
for i in range(len(SubjNamez)):

    NormSpeed = SubjData[SubjNamez[i]+'Norm']['S_AvgSpd']
    NormSpd = np.divide(SubjData[SubjNamez[i]+'Norm']['F_Spd'], NormSpeed)
    
    for t in [0, 1, 2, 3, 4]:
    
        time = SubjData[SubjNamez[i]+Levels[t]]['F_Time']
        TrlSpd = np.divide(SubjData[SubjNamez[i]+Levels[t]]['F_Spd'], NormSpeed)
        # plt.plot(np.std(TrlSpd),'o')
        SpeedVar[i, t] = np.std(TrlSpd)


# speed variability by met metabolic power
X1 = SpeedVar
Y1 = WAvg_F
for i in range(5): # scatter plot all measures
    ax.scatter(X1[:,i], Y1[:,i], 
                c=Colors[i], marker='s', s=sz)
         
# calculate trendlines
z_S = np.polyfit(np.hstack(X1), np.hstack(Y1), 1)
p_S = np.poly1d(z_S)
x = np.linspace(np.min(X1), np.max(X1), 25)
ax.scatter(x,p_S(x),c='k',marker='s', s=sz2, alpha = A)
ax.plot(x,p_S(x),'-k')
# the line equation and R
[c_S, P_S] = stats.pearsonr(np.hstack(X1), np.hstack(Y1))
LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x + ' + str(round(z_S[1],2))
t = ' Fp Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3))
ax.text(0.23, 9, t, fontsize=txt, ha='right', va='top')
ax.set_xlabel('Speed Standard Deviation (m/s)')
ax.set_ylabel('Net Metabolic Power (W/kg)')


# speed variability by CoT
X1 = SpeedVar
Y1 = CoT_Fp_F
for i in range(5): # scatter plot all measures
    ax2.scatter(X1[:,i], Y1[:,i], 
                c=Colors[i], marker='s', s=sz)
         
# calculate trendlines
z_S = np.polyfit(np.hstack(X1), np.hstack(Y1), 1)
p_S = np.poly1d(z_S)
x = np.linspace(np.min(X1), np.max(X1), 25)
ax2.scatter(x,p_S(x),c='k',marker='s', s=sz2, alpha = A)
ax2.plot(x,p_S(x),'-k')
# the line equation and R
[c_S, P_S] = stats.pearsonr(np.hstack(X1), np.hstack(Y1))
LineEq_S = 'y = ' + str(round(z_S[0],2)) + 'x + ' + str(round(z_S[1],2))
t = ' Fp Clamp \n ' + LineEq_S + ' \n R$^2$ = ' + str(round(c_S*c_S,3))
ax2.text(0.23, 5, t, fontsize=txt, ha='right', va='top')
ax2.set_xlabel('Speed Standard Deviation (m/s)')
ax2.set_ylabel('CoT (J/kg/m)')


