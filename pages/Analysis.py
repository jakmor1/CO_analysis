import streamlit as st
import pandas as pd
from io import StringIO
import matplotlib.pyplot as plt
from scipy import signal
import statsmodels.api as sm
import numpy as np
import matplotlib.ticker as ticker
import shutil
import os
import time

#Load data to dataframes for 2 lines
@st.cache_data
def get_data(uploaded_files1,uploaded_files2):
    
    
    control_COs = pd.DataFrame(columns=["LibNum","Chr","Start","Stop","Crossover_sites","Size"])
    mutant_COs = pd.DataFrame(columns=["LibNum","Chr","Start","Stop","Crossover_sites","Size"])
    
    control_size = len(uploaded_files1)
    mutant_size = len(uploaded_files2)
    
    if uploaded_files1:
        
        control_size = len(uploaded_files1)


        control_COs = pd.DataFrame(columns=["LibNum","Chr","Start","Stop","Crossover_sites","Size"])
        for filename in uploaded_files1:
            ret = load_file(filename)
            frames = [control_COs, ret]
            control_COs = pd.concat(frames,ignore_index=True)


    if uploaded_files2:
        mutant_size = len(uploaded_files2)


        mutant_COs = pd.DataFrame(columns=["LibNum","Chr","Start","Stop","Crossover_sites","Size"])
        for filename in uploaded_files2:
            ret = load_file(filename)
            frames = [mutant_COs, ret]
            mutant_COs = pd.concat(frames,ignore_index=True)
        
        
    return control_COs, mutant_COs

#Generate CO frequency plot
@st.cache_data
def CO_frequency_cache(control_COs,mutant_COs,max_chromosome,centromeres,chr_ends,color_line1,color_line2,color_centromere,label_control,label_mutant,legend_outside,tight_layout, slider_range, control_size, mutant_size,file_names):
            time.sleep(2)
            occurrences_control, contr_dict_with_lib = occurrences(control_COs)
            occurrences_mutant, mutant_dict_with_lib = occurrences(mutant_COs)
            
            #Create dict and count number of CO in plants samples
            counts_control = dict()
            for i in occurrences_control:
              counts_control[i] = counts_control.get(i, 0) + 1

            counts_mutant = dict()
            for i in occurrences_mutant:
              counts_mutant[i] = counts_mutant.get(i, 0) + 1
            

            common_keys = set(counts_control.keys()).union(set(counts_mutant.keys()))
            sorted_common_keys = sorted(common_keys)
            
            min_key = int(min(sorted_common_keys))
            max_key = int(max(sorted_common_keys))

            
            counts_control = filter_dict_by_range(counts_control,slider_range[0],slider_range[1])
            counts_mutant = filter_dict_by_range(counts_mutant,slider_range[0],slider_range[1])
            
            lst_removed_lib_control = deleted_libraries(contr_dict_with_lib,slider_range[0],slider_range[1])
            lst_removed_lib_mutant = deleted_libraries(mutant_dict_with_lib,slider_range[0],slider_range[1])
            
            
            #remove chosen lib_num from dataframe
            pd.DataFrame(columns=["LibNum","Chr","Start","Stop","Crossover_sites","Size"])
            
            control_COs = control_COs[~control_COs["LibNum"].isin(lst_removed_lib_control)]
            mutant_COs = mutant_COs[~mutant_COs["LibNum"].isin(lst_removed_lib_mutant)]
            
            #remove deleted files from counting
            control_size = control_size - len(lst_removed_lib_control)
            mutant_size = mutant_size - len(lst_removed_lib_mutant)
            
            #Solving the problem of overlapping numbers on the x-axis by removing adjacent elements 
            common_keys = set(counts_control.keys()).union(set(counts_mutant.keys()))
            sorted_common_keys = sorted(common_keys)
            
            x_ticks = sorted_common_keys 
            x_ticks = remove_adjacent_numbers(sorted_common_keys)
            plt.figure()
            
            #Lower bar on top to make him always visible
            for key in sorted_common_keys:

                control_value = counts_control.get(key, 0)
                mutant_value = counts_mutant.get(key, 0)

                if control_value >= mutant_value:
                    plt.bar(key, control_value, color=color_options[color_line1], align='center', label=label_control if key == sorted_common_keys[0] else '')
                    plt.bar(key, mutant_value, color=color_options[color_line2], align='center',
                            label=label_mutant if key == sorted_common_keys[0] else '')
                else:
                    plt.bar(key, mutant_value, color=color_options[color_line2], align='center', label=label_mutant if key == sorted_common_keys[0] else '')
                    plt.bar(key, control_value, color=color_options[color_line1], align='center',
                            label=label_control if key == sorted_common_keys[0] else '')

            #Counting and visualization of avg values
            suma=0
            for klucz, wartość in counts_control.items():
                suma += klucz * wartość

            plt.axvline(x=suma /(control_size), label = "Avg "+ label_control, color = lighten_color(color_options[color_line1]), linestyle='--' )
            
            suma=0
            for klucz, wartość in counts_mutant.items():
                suma += klucz * wartość
                
            #Plotting
            plt.axvline(x=suma /(mutant_size), label = "Avg " + label_mutant, color = lighten_color(color_options[color_line2]), linestyle='--' )
            plt.xlabel('Crossovers', fontsize=12)
            plt.ylabel('No. of plants', fontsize=12)
            plt.title(label='Crossover frequency', fontsize=20)
            plt.xticks(x_ticks, fontsize = 10)
            
            if legend_outside:
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
            else:
                plt.legend(loc='upper right')
            
            if tight_layout:
                plt.tight_layout()
            plt.savefig('Plot_CO_numbers.png', dpi=1000)
            file_names.append('Plot_CO_numbers.png') 
            st.pyplot(plt.gcf())
            
            plt.close()
            
            #Write removed libraries
            if lst_removed_lib_control:
                st.write("Removed libraries from analisys for "+ label_control+ ":")
                st.write(lst_removed_lib_control) 
            if lst_removed_lib_mutant:
                st.write("Removed libraries from analisys for "+ label_mutant+ ":")
                st.write(lst_removed_lib_mutant)
            
            return control_COs,mutant_COs,lst_removed_lib_control, lst_removed_lib_mutant, control_size, mutant_size,file_names

#Generate CO distribution in chromosomes
@st.cache_data
def CO_chromosome_cache(control_COs, mutant_COs,max_chromosome,centromeres,chr_ends,color_line1,color_line2,color_centromere,window_size,conv_val, advanced_option, file_names,control_size,mutant_size,label_control,label_mutant,legend_outside,tight_layout, max_xsticks=0, max_ysticks=0, insert_size="300 kb"):
    time.sleep(2)
    ma = convolve_value(conv_val) 


    for chromosome in range(1,max_chromosome+1): #for every chromosome
        tmp_control = control_COs[control_COs.Chr == chromosome] #take rows with chosen number of chromosome
        tmp_mutant = mutant_COs[mutant_COs.Chr == chromosome]
        window = list(range(1, chr_ends[chromosome-1], window_size)) #create windows wit chosen step and create last window
        window.append(chr_ends[chromosome-1])

        last_el_in_window = (chr_ends[chromosome-1] - window[len(window)-2]) /window_size #last window normalisation


        window_control = []
        window_mutant = []
        for i in range(len(window)-1): #loop over values in windows and count CO in windows
            COs_in_window = len(tmp_control[(tmp_control.Crossover_sites >= window[i]) & (tmp_control.Crossover_sites <= window[i+1])])
            window_control.append(COs_in_window)
            COs_in_window = len(tmp_mutant[(tmp_mutant.Crossover_sites >= window[i]) & (tmp_mutant.Crossover_sites <= window[i+1])])
            window_mutant.append(COs_in_window)


        window_control[len(window_control)-1] =  window_control[len(window_control)-1]/last_el_in_window
        window_mutant[len(window_mutant) - 1] = window_mutant[len(window_mutant) - 1] / last_el_in_window
        
        #normalisation 
        window_control = [x / control_size for x in window_control]
        window_mutant = [x / mutant_size for x in window_mutant]

        
        #convolving
        if conv_val != 1:
            window_control = modify_list(window_control, conv_val-1)
            window_mutant = modify_list(window_mutant, conv_val-1)

        filt_control = signal.convolve(window_control, ma, mode='valid')
        filt_mutant = signal.convolve(window_mutant, ma, mode='valid')


        plt.plot(window[:-1], filt_control, color=color_options[color_line1], label=label_control)
        plt.plot(window[:-1], filt_mutant, color=color_options[color_line2], label=label_mutant)
        plt.axvline(centromeres[chromosome-1], linestyle='--', color=color_options[color_centromere], label="Centromere")
        
        # Formatting X-axis to wanted order of magnitude 
        plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(format_func))

        if advanced_option and max_xsticks !=0:
            plt.xlim(0,max_xsticks)
            
        if advanced_option and max_ysticks !=0:
            plt.ylim(0,max_ysticks) 

        plt.title("COs distribution along chromosome " + str(chromosome))
        plt.xlabel('Position (Mb)', fontsize=12)
        plt.ylabel('Crossovers per ' + insert_size + ' per F2', fontsize=12)
        if legend_outside:
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        else:
            plt.legend(loc='upper right')
        
        if tight_layout:
            plt.tight_layout()
        plt.savefig('Plot_chromosome_'+str(chromosome)+ '.png', dpi=1000)
        file_names.append('Plot_chromosome_'+str(chromosome) + '.png')
        st.pyplot(plt.gcf())
        plt.close()

    
#Available color options
color_options = {
    "Green": "#008000",
    "Red": "#FF0000",
    "Purple": "#FF00FF",
    "Seagreen3": "#43cd80",
    "Steelblue": "#4682B4",
    "Chocolate": "#d2691e",
    "Burlywood": "#deb887",
    "Mediumpurple": "#9370db",
    "Dimgrey": "#696969",
    "Plum3": "#cd96cd",
    "Lightblue": "#add8e6",
    "Olivedrab3": "#6b8e23",
    "Cyan3": "#00FFFF"

}

def main():
    st.title("Application for analyzing the&nbsp;distribution&nbsp;of crossing-over in&nbsp;Arabidopsis thaliana")
    
    file_names = []
    settings = st.container()
    upload_section = st.container()
    CO_freq = st.container()
    CO_chrom = st.container()
    CO_arms = st.container()
    
    
    lst_removed_lib_control=[]
    lst_removed_lib_mutant =[]
    control_size = 0 
    mutant_size = 0
    centromeres = [15086045,3607929,13587786,3956021,11725024]
    chr_ends = [30427671,19698289,23459830,18585056,26975502]  
    
    color_line1 = "Green"
    color_line2 = "Red"
    color_centromere = "Purple"
    label_control= "Line 1"
    label_mutant= "Line 2"
    
    #General settings
    with settings:
        st.header("General settings:")
        
        #Names, centromeres and length of chromosomes
        col1, col2 = st.columns(2)
        label_control = col1.text_input("Enter name for first line:", value = "Line 1")
        label_mutant = col1.text_input("Enter name for second line:", value = "Line 2")
        centromeres = col1.text_input('Centromeres:', value = "15086045,3607929,13587786,3956021,11725024")
        centromeres = centromeres.split(",")
        centromeres = [eval(i) for i in centromeres]

        chr_ends = col1.text_input('Ends of chromosomes:', value = "30427671,19698289,23459830,18585056,26975502")
        chr_ends = chr_ends.split(",")
        chr_ends = [eval(i) for i in chr_ends]
        
        
        #Color selection
        color_line1 = col2.selectbox("Choose color for first line:", list(color_options.keys()), list(color_options.keys()).index(color_line1))
        color_line2 = col2.selectbox("Choose color for second line:", list(color_options.keys()), list(color_options.keys()).index(color_line2))
        color_centromere = col2.selectbox("Choose color for centromeres:", list(color_options.keys()), list(color_options.keys()).index(color_centromere))


    
    #Upload section
    with upload_section:
        st.header("Upload section:")
        download_button = st.download_button(label='Sample data for analysis', data=open(f'Example_data.zip', 'rb'), file_name='Example_data.zip', mime='application/zip')
        col1, col2 = st.columns(2)
        
        uploaded_files1 = col1.file_uploader("Upload files for first line:", key="control", type=["txt"], accept_multiple_files=True)
        uploaded_files2 = col2.file_uploader("Upload files for second line:",key="mutant", type=["txt"], accept_multiple_files=True)
        
        
        control_COs = pd.DataFrame(columns=["LibNum","Chr","Start","Stop","Crossover_sites","Size"])
        mutant_COs = pd.DataFrame(columns=["LibNum","Chr","Start","Stop","Crossover_sites","Size"])
        
        control_size = len(uploaded_files1)
        mutant_size = len(uploaded_files2)
        

        try:
            control_COs,mutant_COs = get_data(uploaded_files1,uploaded_files2)
        except Exception:
            st.error("Please enter a valid input with correct file format")
            return 0
    
    
    
    if uploaded_files1 and uploaded_files2:

        #CO frequency
        with CO_freq:
            st.header("Analysis of CO frequency:")
            col1, col2 = st.columns(2)
            
                        
            occurrences_control, contr_dict_with_lib = occurrences(control_COs)
            occurrences_mutant, mutant_dict_with_lib = occurrences(mutant_COs)
            
            counts_control = dict()
            for i in occurrences_control:
              counts_control[i] = counts_control.get(i, 0) + 1

            counts_mutant = dict()
            for i in occurrences_mutant:
              counts_mutant[i] = counts_mutant.get(i, 0) + 1
            

            #Set of keys for choosing the range of analyzed CO
            common_keys = set(counts_control.keys()).union(set(counts_mutant.keys()))
            sorted_common_keys = sorted(common_keys)
            
            min_key = int(min(sorted_common_keys))
            max_key = int(max(sorted_common_keys))
            
            
            #Settings for CO frequency
            slider_range = col1.slider("Choose range of crossovers:", value = [min_key,max_key], min_value =min_key, max_value = max_key )
            legend_outside = col2.checkbox('Legend outside of the chart')
            tight_layout = col2.checkbox('Tight layout')
            
            max_chromosome = int(control_COs["Chr"].max())
            
            
            try:
                control_COs,mutant_COs,lst_removed_lib_control, lst_removed_lib_mutant, control_size, mutant_size,file_names = CO_frequency_cache(control_COs,mutant_COs,max_chromosome,centromeres,chr_ends,color_line1,color_line2,color_centromere,label_control,label_mutant,legend_outside,tight_layout, slider_range, control_size, mutant_size,file_names)
            except Exception as err:
                #st.error(f"Unexpected {err=}, {type(err)=}")
                st.error("Wait for the plot to load before you change its variables")
                return 0
    
        
        #Chromosome plots
        #Decreasing value - eliminated libraries in settings
        control_size = control_size - len(lst_removed_lib_control)
        mutant_size = mutant_size - len(lst_removed_lib_mutant)
            
        with CO_chrom:
            st.header("Analysis of CO of individual chromosomes:")
            col1, col2 = st.columns(2)
            
            #Settings for chromosomes plots
            conv_val = col1.slider("Choose the degree of convolution:", min_value=1,max_value=31, value = 7, step = 2)
            window_size = col1.number_input('Insert a window size in bp:', value = 300000, min_value = 1, max_value = max(chr_ends))
            advanced_option = col2.checkbox('Show advanced options for plotting')
            
            if advanced_option:
                max_xsticks = col2.number_input("Insert max value for X-axis in Mb:", min_value = 0.0, value =0.0)
                max_xsticks = max_xsticks * 1000000 #In Mb
                max_ysticks = col2.number_input("Insert max value for Y-axis:", min_value = 0.0, value =0.0)
                insert_size = col2.text_input("Enter window size for title in Y-axis:", value = "300 kb")
            
                try:
                    CO_chromosome_cache(control_COs, mutant_COs,max_chromosome,centromeres,chr_ends,color_line1,color_line2,color_centromere,window_size,conv_val, advanced_option, file_names,control_size,mutant_size,label_control,label_mutant,legend_outside,tight_layout, max_xsticks, max_ysticks, insert_size)
                except IndexError as err:
                    st.error("Check if the number of analyzed chromosomes matches the number of provided centromeres or chromosome ends." )
                    return 0
                
                except Exception as err:
                    st.error(f"Unexpected {err=}, {type(err)=}")
                    st.error("Wait for the plots to load before you change its variables ver1" )
                    return 0
            else:
                try:
                    CO_chromosome_cache(control_COs, mutant_COs,max_chromosome,centromeres,chr_ends,color_line1,color_line2,color_centromere,window_size,conv_val, advanced_option,file_names,control_size,mutant_size,label_control,label_mutant,legend_outside,tight_layout)
                except IndexError as err:
                    st.error("Check if the number of analyzed chromosomes matches the number of provided centromeres or chromosome ends." )
                    return 0
                except Exception as err:
                    st.error(f"Unexpected {err=}, {type(err)=}")
                    st.error("Wait for the plots to load before you change its variables ver2" )
                    return 0



        #TEL-CEN plotting
        with CO_arms:
            st.header("Analysis of CO along chromosome arms:")
            col1, col2 = st.columns(2)
            
            control_all_prob = pd.DataFrame()
            mutant_all_prob= pd.DataFrame()


            wins = [x/1000 for x in range(0, 1001)]

            wins_control = []
            wins_mutant = []
            #Prepare data for plotting
            for chromosome in range(1,max_chromosome+1):
                tmp_all_prob = pd.DataFrame()

                tmp_all_prob = tel_cen_data(control_COs,chromosome,centromeres,chr_ends)
                control_all_prob= pd.concat([control_all_prob, tmp_all_prob],ignore_index=True)

                tmp_all_prob = tel_cen_data(mutant_COs, chromosome,centromeres,chr_ends)
                mutant_all_prob= pd.concat([mutant_all_prob, tmp_all_prob],ignore_index=True)

            control_all_prob = control_all_prob.rename(columns={control_all_prob.columns[0]: 'Stop'})
            mutant_all_prob = mutant_all_prob.rename(columns={mutant_all_prob.columns[0]: 'Stop'})


            for i in range(len(wins)-1):
                tmp_win = len(control_all_prob[(control_all_prob.Stop>=wins[i]) & (control_all_prob.Stop<wins[i+1])])
                wins_control.append(tmp_win)

                tmp_win = len(mutant_all_prob[(mutant_all_prob.Stop>=wins[i]) & (mutant_all_prob.Stop<wins[i+1])])
                wins_mutant.append(tmp_win)

            #Normalisation
            control_norm = [x / control_size for x in wins_control]
            mutant_norm = [x / mutant_size for x in wins_mutant]


            xlim = [0, 1]      
            wins2 = wins[1:]

            fig, ax = plt.subplots()
            ax.set_xlim(xlim)   
            ax.set_title("TEL-CEN plotting")
            
            #plot data and smooth spline
            x = np.array(wins2)
            y1 = np.array(control_norm)
            y2 = np.array(mutant_norm)

            frac_number = col1.slider("Choose fraction size (degree of smoothing) in LOWESS:", min_value=0.00,max_value=0.5, value = 0.08, step = 0.005)
            
            lowess_model = sm.nonparametric.lowess(y1, x, frac = frac_number) 


            plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1], ["TEL", 0.2, 0.4, 0.6, 0.8, "CEN"])
            max_ysticks = col2.number_input("Insert max value for Y-axis:",key ="TEL-CEN", min_value = 0.0, value =0.0,format="%.3f")

            maxim = max(lowess_model[:, 1])
            
            # Plotting
            plt.plot(lowess_model[:, 0], lowess_model[:, 1],  label=label_control, color=color_options[color_line1], zorder=3)
            
            ax.axhline(y=np.mean(lowess_model[:, 1]), color=lighten_color(color_options[color_line1]), linestyle='--', label="Avg " + label_control)
            
            lowess_model = sm.nonparametric.lowess(y2, x, frac = frac_number)  

            plt.plot(lowess_model[:, 0], lowess_model[:, 1], label=label_mutant, color=color_options[color_line2], zorder=3)
            
            if(max(lowess_model[:, 1]) > maxim):
                maxim = max(lowess_model[:, 1])
            
            plt.ylim(0, maxim *1.1)
            
            if max_ysticks !=0:
                plt.ylim(0,max_ysticks) 
            
            ax.axhline(y=np.mean(lowess_model[:, 1]), color=lighten_color(color_options[color_line2]), linestyle='--', label="Avg " + label_mutant)
            
            if legend_outside:
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
            else:
                plt.legend(loc='upper right')
                  
            
            plt.title("COs distribution along chromosome arms")
            plt.xlabel('Position', fontsize=12)
            plt.ylabel('Crossovers', fontsize=12)
            if tight_layout:
                plt.tight_layout()
            plt.savefig('Plot_TEL-CEN.png', dpi=1000)
            file_names.append('Plot_TEL-CEN.png')
            st.pyplot(plt.gcf())
            plt.close()
            
            
            
            #Save Button
            # Create a temporary folder to save files
            temp_folder = 'temp_images'
            os.makedirs(temp_folder, exist_ok=True)
            
            # Copy files to temporary folder
            for file_name in file_names:
                shutil.copy(file_name, os.path.join(temp_folder, file_name))
            
            # Create ZIP
            shutil.make_archive(temp_folder, 'zip', temp_folder)
            
            # PDownload ZIP
            st.download_button(label='Download generated plots', data=open(f'{temp_folder}.zip', 'rb'), file_name='CO_plots.zip', mime='application/zip')
            
            # Remove temporary folder and ZIP
            shutil.rmtree(temp_folder)
            os.remove(f'{temp_folder}.zip')

 

def remove_adjacent_numbers(lst): #remove adjacent numbers in the array
    result = []
    for num in lst:
        if not result or abs(num - result[-1]) > 2: #by at least 2 e.g. [0,1,2,3,4] = [0,3]
            result.append(num)
    return result


def load_file(file_name): #Load data from files

    str = file_name.name
    str = str.split(".")
    str = str[1]
    

    stringio = StringIO(file_name.getvalue().decode("utf-8"))

    my_file=stringio.read()

    
    lines = my_file.split("\n")
    lst = []
    for line in lines:

        lst.append(line.split())

    lst = [[int(x) if x.isdigit() else x for x in subarray] for subarray in lst]  #convert string to numbers

    df_datafile = pd.DataFrame(data=lst, columns=["Col0", "Col1", "Col2", "Col3", "Col4"])
    
    max_chromosome = int(df_datafile["Col1"].max())

    dane = []
    for chromosome in range(1, max_chromosome+1):
        result = df_datafile[df_datafile["Col1"] == chromosome]
        if (len(result.index)> 1): #more than 1 genotype in chromosome

            starts = result["Col3"].tolist()[:-1]  # without the last locus termination point
            ends = result["Col2"].tolist()[1:]  # without the first locus termination point
            num = int(str)

            for i in range(len(starts)):
                COs = starts[i] + ((ends[i] - starts[i]) // 2)
                size = ends[i] - starts[i]
                tmp_lst = [num, chromosome, starts[i], ends[i], COs, size]
                dane.append(tmp_lst)

    df_perfile = pd.DataFrame(data=dane, columns=["LibNum", "Chr", "Start", "Stop", "Crossover_sites", "Size"])

    return df_perfile

def occurrences(df): #count crossing-over occurences for every unique lib(column "LibNum")
    lst = df["LibNum"].tolist()
    lst = list(set(lst))

    dict = {}

    occurrences = []
    for unique_num in lst:
        count = df["LibNum"].value_counts()[unique_num]
        occurrences.append(count)
        dict[unique_num] = count

    return occurrences, dict
    
#create list to convolve
def convolve_value(x):
    if x <= 0:
        return []

    return [1 / x] * x
    
#Prepare data to TEL-CEN analysis
def tel_cen_data(control_COs,chromosome,centromeres,chr_ends):
    tmp_control = control_COs[control_COs.Chr == chromosome]
    tmp_control = tmp_control[["Stop"]]

    control_left = tmp_control[tmp_control.Stop < centromeres[chromosome-1]] 
    control_left_prop = control_left["Stop"].div(centromeres[chromosome-1])


    control_right = tmp_control[tmp_control.Stop > centromeres[chromosome-1]] 
    control_right["Stop"] = control_right["Stop"] - centromeres[chromosome-1]

    control_right_prop = pd.DataFrame(columns = ["Stop"])
    control_right_prop["Stop"] = control_right["Stop"].div((chr_ends[chromosome-1] - centromeres[chromosome-1]))  


    control_right_prop["Stop"] = control_right_prop["Stop"].div(-1)
    control_right_prop["Stop"] = control_right_prop["Stop"] + 1  
    control_right_prop = control_right_prop.Stop[::-1]  # reversed
    control_right_prop.reset_index(inplace=True, drop=True)


    control_all_prob = pd.concat([control_left_prop, control_right_prop],ignore_index=True)
    return control_all_prob

#Normalisation of convolving edges
def modify_list(arr,number):
    number = number//2
    first_x = arr[:number]  # Copy the first 4 elements
    last_x = arr[-number:]  # Copy the last 4 elements

    #Inserting
    arr = first_x + arr + last_x

    return arr
    
#Filter by range of crossovers
def filter_dict_by_range(dictionary, min_range, max_range):
    keys_to_remove = []
    
    for key in dictionary.keys():
        if key < min_range or key > max_range:
            keys_to_remove.append(key)
    
    for key in keys_to_remove: 
        del dictionary[key]
    
    return dictionary
    
#Remember deleted libraries
def deleted_libraries(dictionary, min_range, max_range):
    removed_lib = []
    
    for key, value in dictionary.items():
        if value < min_range or value > max_range:
            removed_lib.append(key)
    
    return removed_lib

#Pair base to Mb
def format_func(value, tick_number):
    return int(value / 1000000)
    
hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

#Lighten color
def lighten_color(hex_color, factor=0.5):

    

    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (1, 3, 5))


    new_r = int(r + (255 - r) * factor)
    new_g = int(g + (255 - g) * factor)
    new_b = int(b + (255 - b) * factor)


    new_hex_color = "#{:02X}{:02X}{:02X}".format(new_r, new_g, new_b)

    return new_hex_color

if __name__ == "__main__":
    main()
