import numpy as np 
import re 
import subprocess 
from collections import Counter
import os 
import os.path

pw_file='01-relax-pmwo'

pw_file_in=pw_file+".in"
pw_file_out=pw_file+'.out'

cwd=os.getcwd()
path=cwd + '/' + pw_file_in

print(path)
if os.path.isfile(path):
   print("*******Input file name********", pw_file_in)
else:
    print("*******Input file does not exit********")
    pw_file_in=input(" Enter input filename\n")



def vasp_output(coord_sys, natms_names, natms, ntype_atms, latt_vec, pos):
    fout_vasp = open("struc.vasp", 'w')

    # Cartesian to Direct
    if coord_sys=="angstrom":
       pos= np.matmul(np.copy(pos), np.linalg.inv(latt_vec))
    

    fout_vasp.writelines(natms_names)
    fout_vasp.write("\n")
    fout_vasp.write("1.0" + "\n")
    for i in range(3):
        fout_vasp.write('{0:12.8f} {1:12.8f} {2:12.8f}'.format(latt_vec[i][0], latt_vec[i][1], latt_vec[i][2]) + "\n")

    nspecies=int(len(natms_names))

    for i in range(0, nspecies):
        fout_vasp.write('{:s}'.format(natms_names[i])+ '  ')
    fout_vasp.write("\n")
    for i in range(0, nspecies):
        fout_vasp.write('{:d}'.format(ntype_atms[i]) + '  ')
    fout_vasp.write("\n")
    if coord_sys=='crystal':
       fout_vasp.write("Direct" + "\n")
    elif coord_sys=='angstrom':
       fout_vasp.write("Cartesian" + "\n")

    for i in range(0, natms):
        x=pos[i][0]
        y=pos[i][1]
        z=pos[i][2]
        fout_vasp.write('{:12.8f} {:12.8f} {:12.8f}'.format(x, y, z) + "\n")

def qe_imput(fin, natms):
    #-----------------------qe-input-------------------------------------
    command='grep -A4 "CELL_PARAMETERS angstrom" ' + fin +  ' |' + " awk '{print $1, $2, $3}' | sed -e '1,1d' "
    p1=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    output, error=p1.communicate()
    p1.kill()

    out1= output.split('\n')

    cell=np.zeros((3,3))


    for i in range(0, 3):
        x=' '.join(out1[i].split()).split() # remove epmpty space and put into a list
        #print(x)
        cell[i,:]=np.array(x,dtype=np.float64)
    del x
    command='grep -A' + str(natms) + ' "ATOMIC_POSITIONS" ' + fin + ' |' + " awk '{print $1}'  | sed -e '1,1d' "
    p1=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    output, error=p1.communicate()
    p1.kill()

    out1= output.split('\n')
    #print(out1[0:natms])
    no_names=out1[0:natms]
    #ntype_names=list(set(natms_names))
    natms_names=list(Counter(no_names).keys()) # equals to list(set(words))
    ntype_atms=list(Counter(no_names).values()) # counts the elements' frequency
    
    command='grep -A' + str(natms) + ' "ATOMIC_POSITIONS" ' + fin + ' |' + " awk '{print $2, $3, $4}'  | sed -e '1,1d' "
   # print(command)
    p1=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    output, error=p1.communicate()
    p1.kill()

    out1= output.split('\n')
    pos_crystal=np.zeros((natms,3), dtype=float)
     
    for i in range(0, natms):
        x=' '.join(out1[i].split()).split() # remove epmpty space and put into a list
        #print(x)
        #temp_forces=np.array(x[1:len(x)],dtype=np.float64)
        pos_crystal[i,:]=np.array(x,dtype=np.float64)

    del out1
    del x
    return (natms_names, ntype_atms, cell, pos_crystal)


def qe_imput_cell(fin, natms):
    #-----------------------qe-input-------------------------------------
    command='grep -A3 "CELL_PARAMETERS" ' + fin +  ' |' + " awk '{print $1, $2, $3}' | tail -3 " # sed -e '1,1d' "
    p1=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    output, error=p1.communicate()
    p1.kill()

    out1= output.split('\n')

    cell=np.zeros((3,3))

    print(out1)
    for i in range(0, 3):
        x=' '.join(out1[i].split()).split() # remove epmpty space and put into a list
        #print(x)
        cell[i,:]=np.array(x,dtype=np.float64)
    del x
    del out1
    print(cell)
    return cell
#natms=40

def qe_imput_pos(fin, natms):
    command='grep -A' + str(natms) + ' "ATOMIC_POSITIONS" ' + fin + ' |' + " awk '{print $1}'  | tail -" + str(natms)  #sed -e '1,1d' "
    print(command)

    p1=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    output, error=p1.communicate()
    p1.kill()

    out1= output.split('\n')
    #print(out1[0:natms])
    no_names=out1[0:natms]
    #ntype_names=list(set(natms_names))
    natms_names=list(Counter(no_names).keys()) # equals to list(set(words))
    ntype_atms=list(Counter(no_names).values()) # counts the elements' frequency

    command='grep -A' + str(natms) + ' "ATOMIC_POSITIONS" ' + fin + ' |' + " awk '{print $2, $3, $4}'  | tail -" + str(natms) #sed -e '1,1d' "
   # print(command)
    p1=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
    output, error=p1.communicate()
    p1.kill()

    out1= output.split('\n')
    pos_crystal=np.zeros((natms,3), dtype=float)
#    print(out1)

    for i in range(0, natms):
        x=' '.join(out1[i].split()).split() # remove epmpty space and put into a list
        #print(x)
        #temp_forces=np.array(x[1:len(x)],dtype=np.float64)
        pos_crystal[i,:]=np.array(x,dtype=np.float64)

    del out1
    del x
    return (natms_names, ntype_atms,pos_crystal)

with open(pw_file_in, 'r') as inpw:
    for line in inpw:
        line = line.rstrip()
        # search case insensative 
        if re.search('calculation', line, flags=re.IGNORECASE):
            # string after equals to sign
            s=line.split("=")[1].strip() # 
            # Remove specific characters from a string
            cal_type = re.sub("[!@#$']", '', s) # remove character any of [!@#$] from word
            #print(line)
        if re.search('nat', line, flags=re.IGNORECASE):
            # string after equals to sign
            s=line.split("=")[1].strip() # 
            # Remove specific characters from a string
            natms = int(re.sub("[!@#$']", '', s)) # remove character any of [!@#$] from word
            print(natms)

        if re.search('ntyp', line, flags=re.IGNORECASE):
            # string after equals to sign
            s=line.split("=")[1].strip() # 
            # Remove specific characters from a string
            ntypes = int(re.sub("[!@#$']", '', s)) # remove character any of [!@#$] from word
            print(ntypes)

        if re.search('ATOMIC_POSITIONS', line, flags=re.IGNORECASE):
            # string after equals to sign
            s=line.split(" ")[1].strip() # 
            # Remove specific characters from a string
            coord_sys =re.sub("[(!@#$')]", '', s) # remove character any of [!@#$] from word
            print(coord_sys)


def cehck_conv(pw_file_out):
    with open(pw_file_out, 'r') as inpw:
        for line in inpw:
            line = line.rstrip()
        if re.search('Begin', line):
            print("structure converged")
        else:
            print("structure do not converge")

    return

    #tsteps=[int(i) for i in x.split()]

#latt_vec=qe_imput_cell(pw_file_in, natms)
#natms_names, ntype_atms, pos=qe_imput_pos(pw_file_in, natms)
#print(pos)


path=cwd + '/' + pw_file_out
#print(path)
if os.path.isfile(path):
   print("*******output file name********", pw_file_out)
else:
    print("*******output file does not exit********")
    pw_file_out=input(" Enter output filename\n")

if cal_type=='relax':
    path=cwd + '/' + pw_file_out
    latt_vec=qe_imput_cell(pw_file_in, natms)
    natms_names, ntype_atms, pos=qe_imput_pos(pw_file_out, natms)
    cehck_conv(pw_file_out)

elif cal_type=='vc-relax':
    print('vc-relax')
    latt_vec=qe_imput_cell(pw_file_out, natms)
    natms_names, ntype_atms, pos=qe_imput_pos(pw_file_out, natms)
    cehck_conv(pw_file_out)
else:
    natms_names, ntype_atms, latt_vec, pos=qe_imput(pw_file_in, natms)

vasp_output(coord_sys, natms_names, natms, ntype_atms, latt_vec, pos)

