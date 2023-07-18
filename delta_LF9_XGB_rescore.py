import pandas as pd

df=pd.read_csv('poses_data.csv')
IScore=pd.DataFrame()
poses_id,lig,o_index,Lin_F9_Score=[],[],[],[]
for i in list(df.poses_id.to_list()):
    print(i)
    lig_file='/vast/zf2012/diffdock/Mtb_folder_combined/poses/%s.pdb'%i
    f1=open(lig_file,'r')
    for line in f1.readlines():
        if line.startswith('REMARK minimizedAffinity'):
            score=line.split(' ')[2].split('\n')[0]
            Lin_F9_Score.append(score)
            poses_id.append(i)
            lig.append(i.split('_')[0])
            o_index.append(i.split('_')[1])
IScore['poses_id']=poses_id
IScore['lig']=lig
IScore['lig']=IScore['lig'].astype(str)
IScore['o_index']=o_index
IScore['Lin_F9_Score']=Lin_F9_Score

df2=pd.read_csv('../combined_smiles_Class.csv')
df2['lig']=df2['Index'].astype(str)
left=IScore.set_index('lig')
right=df2.set_index('lig')
df_LF9=left.join(right)
df_LF9.to_csv('poses_data_Lin_F9_Score.csv')

import sys,os
import numpy as np
from scipy.spatial.distance import cdist

def get_lig_coord(filename):
    '''
    Get all heavy atom coordinates and return a coordinate matrix

    '''
    f1 = open(filename, 'r')

    List = []
    for line in f1.readlines():
        if line.startswith(('HETATM', 'ATOM')) and line[77]!='H':
            coord_x, coord_y, coord_z = line[27:38], line[38:46], line[46:54]
            result = [float(coord_x), float(coord_y), float(coord_z)]
            List.append(result)
    f1.close()
    List = np.asarray(List)
    return List

def get_beta_info(filename):
    '''
    Get beta atom coordinates and the corresponding beta atom scores

    '''
    f1 = open(filename, 'r')

    List = []
    Score = []
    for line in f1.readlines():
        if line.startswith(('HETATM', 'ATOM')) and line[77]!='H':
            coord_x, coord_y, coord_z = line[27:38], line[38:46], line[46:54]
            result = [float(coord_x), float(coord_y), float(coord_z)]
            List.append(result)
            score = float(line[61:66])
            print(score)
            Score.append(score)
    f1.close()
    List = np.asarray(List)
    Score = np.asarray(Score)
    return List, Score

def calc_betaScore_and_ligCover(lig_file, beta_file):
    lig_coords = get_lig_coord(lig_file)
    beta_coords, beta_scores = get_beta_info(beta_file)
    result1 = cdist(beta_coords, lig_coords)
    score = np.sum(beta_scores[np.min(result1, axis=1)<=1.6])
    result2 = cdist(lig_coords, beta_coords)
    lig_cover_coords = lig_coords[np.min(result2, axis=1)<=1.6]
    return round(score,3),round(len(lig_cover_coords)/len(lig_coords),3)

lig_Cover,betaScore=[],[]
for i in list(df_LF9.poses_id.to_list()):
    beta_file='/home/zf2012/Mtb_DnaK_AF2_rank0/ranked_0_protein_betaAtoms.pdb'
    lig_file='/vast/zf2012/diffdock/Mtb_folder_combined/poses/%s.pdb'%i
    a1,a2 = calc_betaScore_and_ligCover(lig_file,beta_file)
    betaScore.append(a1)
    lig_Cover.append(a2)

df_LF9['betaScore']=betaScore
df_LF9['lig_Cover']=lig_Cover
df_LF9.to_csv('poses_data_Lin_F9_Score_betaScore_ligCover.csv')


import sys
import rdkit
import random
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

Lin_F9 = round(df_LF9['Lin_F9_Score'].astype(float)*(-0.73349),3)
df_LF9['Lin_F9']=Lin_F9
LigandDescriptors=['HeavyAtomMolWt','NumValenceElectrons','FpDensityMorgan1',
                     'FpDensityMorgan2','FpDensityMorgan3','LabuteASA',
                     'TPSA','NHOHCount','MolLogP','MolMR']
DescCalc = MolecularDescriptorCalculator(LigandDescriptors)
def GetRDKitDescriptors(smile):
    mol = Chem.MolFromSmiles(smile)
    mol.UpdatePropertyCache(strict = False)
    Chem.GetSymmSSSR(mol)
    return DescCalc.CalcDescriptors(mol)

Features = []
print("df_LF9_len:",len(df_LF9))

for i in list(df_LF9['poses_id'].tolist()):
    smile = df_LF9.loc[df_LF9['poses_id']==i, 'smiles'].values[0]
    Features.append(GetRDKitDescriptors(smile))
    #print(i)

ss = pd.DataFrame(Features, columns=LigandDescriptors)
print("ss_len:",len(ss))
print(df_LF9)
print(ss)
df_LF9=df_LF9.reset_index()
df_rdkit_descriptors=pd.concat([df_LF9,ss],axis=1, ignore_index=True)
df_rdkit_descriptors.columns = list(df_LF9.columns) + list(LigandDescriptors)
df_rdkit_descriptors.to_csv('df_rdkit_descriptors.csv')

List = ['metal%d'%x for x in range(2,8)] + ['Nbw','Epw','Elw']
def Update(df,List):
    for i in List:
        df[i] = pd.Series([0]*len(df), index=df.index)
    return df
df_update = Update(df_rdkit_descriptors,List)
df_update.to_csv('poses_data_update.csv')

Vina=['vina%d'%i for i in range(1,49)]
f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
SASA = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]

vina_file=pd.read_csv('poses_data_feature.csv')
sasa_file=pd.read_csv('poses_data_sasa.csv')

left=df_update.set_index('poses_id')
right=vina_file.set_index('poses_id')
common1=left.join(right)

third=sasa_file.set_index('poses_id')
input_data=common1.join(third)
input_data.to_csv('poses_data_all.csv')

import xgboost as xgb
from sklearn.linear_model import LinearRegression
from scipy import stats
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
import matplotlib.pyplot as plt
import pickle


print(xgb.__version__)

f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
SASA = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]
SASA = [x for x in SASA if x not in ['P2dl.HA','P2dp.HA']]

vina1 = ['vina%d'%x for x in range(2,9)]
vina2 = ['vina%d'%x for x in range(10,17)]
vina3 = ['vina%d'%x for x in range(18,25)]
vina4 = ['vina%d'%x for x in range(26,31)]
vina5 = ['vina%d'%x for x in range(32,37)]
vina6 = ['vina%d'%x for x in range(38,49)]

vina=vina1+vina2+vina3+vina4+vina5+vina6
metal = ['metal%d'%x for x in range(2,8)]
BW = ["Nbw","Epw","Elw"]

columns = vina+['lig_Cover','betaScore','LE']+SASA+metal+LigandDescriptors+BW
input_data['LE'] = input_data['Lin_F9']/input_data['vina45']
len(columns)

X = np.c_[input_data[columns]]
X = X.astype(np.float64)
y_fix = np.r_[input_data['Lin_F9']]
y_predict_ = []
for i in range(1,11):
    xgb_model = pickle.load(open("/home/zf2012/delta_LinF9_XGB/XGB_model/mod_12_%d.pickle.dat"%i,"rb"))
    y_i_predict = xgb_model.predict(X, ntree_limit=xgb_model.best_ntree_limit)
    y_predict_.append(y_i_predict)

y_predict = np.average(y_predict_, axis=0)
print('y_predict:',y_predict)
input_data['XGB'] = pd.Series(y_predict+y_fix, index=input_data.index)
print('XGB:',input_data['XGB'])
input_data['XGB_score'] = round(input_data['XGB']/(-0.73349),3)
print('XGB_score:',input_data['XGB_score'])
input_data['XGB_score_LE'] = round(input_data['XGB_score']/input_data['vina45'],3)
print('XGB_score_LE:',input_data['XGB_score_LE'])
input_data['Lin_F9_Score']=input_data['Lin_F9_Score'].astype(float)
print('Lin_F9_Score:',input_data['Lin_F9_Score'])
input_data['LinF9_score_LE'] = round(input_data['Lin_F9_Score']/input_data['vina45'],3)
print('LF9_score_LE:',input_data['LinF9_score_LE'])
input_data.to_csv('output_dataset.csv')

output_5pose = input_data[['lig','smiles','o_index','betaScore','lig_Cover','vina45','Lin_F9_Score','LE','XGB','XGB_score','XGB_score_LE','LinF9_score_LE']].sort_values(by=['lig','XGB'], ascending=[True,False])
output_top1pose=output_5pose.drop_duplicates('lig')
output_top1pose.to_csv('output_top1pose.csv')








