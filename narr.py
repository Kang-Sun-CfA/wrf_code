# -*- coding: utf-8 -*-
"""
Created on Wed May 29 20:07:16 2019

@author: kangsun

python class handling NARR
"""
import datetime
import numpy as np
import os
import sys
import logging
import scipy.io as sio
from scipy.interpolate import RegularGridInterpolator
import requests
from calendar import monthrange

def F_show_grib_variables(fn):
    """
    print out all variables in a grib file
    2019/05/03
    """
    import Nio
    f = Nio.open_file(fn,format='grib')
    varnames = f.variables.keys()
    for var in varnames:
        print(var)
    
def F_open_narr_grib(fn,narr_index_list,narr_name_list,if_save_coordinates=False):
    """
    open narr grib file, PyNIO works
    fn:
        full path of grib file name
    narr_index_list:
        variable names in grib files, try F_show_grib_variables(fn) to see all 
        possible variable names used in a grib file, fn
    narr_name_list:
        what you want to call those variables, have to be the same size as narr_index_list
    if_save_coordinate:
        True or False. It is assumed that coordinates the same for all narr files
    extra packages:
        PyNIO, libiconv
    updated on 2019/05/01, tried pygrib with limited success but cannot handle flx files
    overhaul on 2019/05/02 to switch to PyNIO
    """
    import Nio
    f = Nio.open_file(fn,format='grib')
    narr = {}
    for ivar in range(len(narr_index_list)):
        narr_index = narr_index_list[ivar]
        narr_name = narr_name_list[ivar]
        value = f.variables[narr_index][:]
        if len(value.shape) == 3:# move vertical dimension to the end, after lat/lon
            value = np.transpose(value,(1,2,0))
        narr[narr_name] = np.asfortranarray(value)
        if ('lat' not in narr.keys()) and if_save_coordinates:
            narr['lat'] = np.asfortranarray(f.variables['gridlat_221'][:])
            narr['lon'] = np.asfortranarray(f.variables['gridlon_221'][:])
    return narr

def check_file_status(filepath, filesize):
    sys.stdout.write('\r')
    sys.stdout.flush()
    size = int(os.stat(filepath).st_size)
    percent_complete = (size/filesize)*100
    sys.stdout.write('%.3f %s' % (percent_complete, '% Completed'))
    sys.stdout.flush()

def F_subset_narr_variable(v,startx_index_1based,endx_index_1based,\
                 starty_index_1based,endy_index_1based):
    """
    subset narr variable
    """
    if v.ndim == 1:
        if len(v) == 349:
            v = v[startx_index_1based-1:endx_index_1based]
        elif len(v) == 277:
            v = v[starty_index_1based-1:endy_index_1based]
        else:
            logging.error('this should not happen!')
    elif v.ndim == 2:
        v = v[starty_index_1based-1:endy_index_1based,\
                    startx_index_1based-1:endx_index_1based]
    elif v.ndim == 3:
        v = v[starty_index_1based-1:endy_index_1based,\
                    startx_index_1based-1:endx_index_1based,...]
    return v

class narr(object):
    
    def __init__(self,narr_dir,narr_constants_path='narr_constants.mat',\
                 startx_index_1based=1,endx_index_1based=349,\
                 starty_index_1based=1,endy_index_1based=277):
        """
        initiate the narr object
        """
        self.logger = logging.getLogger(__name__)
        self.logger.info('creating an instance of narr logger')
        if not os.path.exists(narr_dir):
            self.logger.warning('narr_dir '+narr_dir+' does not exist! creating one...')
            os.makedirs(narr_dir)
        self.narr_dir = narr_dir;
        self.startx_index_1based = startx_index_1based
        self.endx_index_1based = endx_index_1based
        self.starty_index_1based = starty_index_1based
        self.endy_index_1based = endy_index_1based
        constants = sio.loadmat(narr_constants_path)
        x0 = constants['x'].flatten()
        y0 = constants['y'].flatten()
        hgt0 = constants['hgt'].squeeze()
        lon0 = constants['lon'].squeeze()
        lat0 = constants['lat'].squeeze()
        x = x0[startx_index_1based-1:endx_index_1based]
        y = y0[starty_index_1based-1:endy_index_1based]
        hgt = F_subset_narr_variable(hgt0,startx_index_1based,endx_index_1based,\
                 starty_index_1based,endy_index_1based)
        lon = F_subset_narr_variable(lon0,startx_index_1based,endx_index_1based,\
                 starty_index_1based,endy_index_1based)
        lat = F_subset_narr_variable(lat0,startx_index_1based,endx_index_1based,\
                 starty_index_1based,endy_index_1based)
        self.lat = lat
        self.lon = lon
        self.lon0 = lon0
        self.lat0 = lat0
        self.x = x
        self.y = y
        self.x0 = x0
        self.y0 = y0
        self.hgt = hgt
        self.hgt0 = hgt0
        self.nx = len(self.x)
        self.ny = len(self.y)
        
        step_hour = 3
        daily_start_time = datetime.time(hour=0,minute=0)
        
        self.step_hour = step_hour
        self.daily_start_time = daily_start_time
        self.nlayer = 30
        self.nlayer0 = 30
        self.nx0 = 349
        self.ny0 = 277
    
    def F_set_time_bound(self,start_year=1995,start_month=1,start_day=1,\
                 start_hour=0,start_minute=0,start_second=0,\
                 end_year=2025,end_month=12,end_day=31,\
                 end_hour=0,end_minute=0,end_second=0):
        """
        reset start and end time.
        also create narr time stamps covering the time bounds
        """
        self.start_python_datetime = datetime.datetime(start_year,start_month,start_day,\
                                                  start_hour,start_minute,start_second)
        self.end_python_datetime = datetime.datetime(end_year,end_month,end_day,\
                                                end_hour,end_minute,end_second)
        step_hour = self.step_hour
        daily_start_time = self.daily_start_time
        # extend the start/end datetime to the closest step_hour intervals
        t_array0 = datetime.datetime.combine(datetime.date(start_year,start_month,start_day),\
        daily_start_time)-datetime.timedelta(hours=step_hour)
        t_array = np.array([t_array0+datetime.timedelta(hours=int(step_hour)*i) for i in range(int(24/step_hour+2))])
        tn_array = np.array([(self.start_python_datetime-dt).total_seconds() for dt in t_array])
        narr_start_datetime = t_array[tn_array >= 0.][-1]
        
        t_array0 = datetime.datetime.combine(datetime.date(end_year,end_month,end_day),\
        daily_start_time)-datetime.timedelta(hours=step_hour)
        t_array = np.array([t_array0+datetime.timedelta(hours=int(step_hour)*i) for i in range(int(24/step_hour+2))])
        tn_array = np.array([(self.end_python_datetime-dt).total_seconds() for dt in t_array])
        narr_end_datetime = t_array[tn_array <= 0.][0]
        
        nstep = (narr_end_datetime-narr_start_datetime).total_seconds()/3600/step_hour+1
        
        self.narr_start_datetime = narr_start_datetime
        self.narr_end_datetime = narr_end_datetime
        self.nstep = int(nstep)
        self.logger.info('specified time from '+\
                         self.start_python_datetime.strftime('%Y-%m-%dT%H:%M:%SZ')+
                         ' to '+self.end_python_datetime.strftime('%Y-%m-%dT%H:%M:%SZ'))
        self.logger.info('extended time from '+\
                         self.narr_start_datetime.strftime('%Y-%m-%dT%H:%M:%SZ')+
                         ' to '+self.narr_end_datetime.strftime('%Y-%m-%dT%H:%M:%SZ'))
        self.logger.info('there will be %d'%nstep+' narr time steps')
    
    def F_generate_filelist(self,file_collection_names):
        """
        generate a list of tar files to be downloaded from ncar rda website
        file_collection_names:
            a list of narr files, chosen from ['3D','sfc','flx','clm','pbl']
            y/m/d should be just numbers. the end date is excluded
        updated on 2019/05/29
        """
        self.file_collection_names = file_collection_names
        collection_list = ['3D','sfc','flx','clm','pbl']
        date_start_list = [np.arange(1,28,3),np.array([1,10,20]),np.array([1,9,17,25]),\
                           np.array([1]),np.array([1,10,20])]
        start_date = self.narr_start_datetime.date()
        end_date = self.narr_end_datetime.date()
        
        days = (end_date-start_date).days
        DATES = [start_date + datetime.timedelta(days=d) for d in range(days)]
        fn_list = []
        for DATE in DATES:
            ystr = "%04d"%DATE.year
            mstr = "%02d"%DATE.month
            for file_collection_name in file_collection_names:
                date_start = date_start_list[collection_list.index(file_collection_name)]
                if DATE.day in date_start[0:-1]:
                    d1str = "%02d"%DATE.day
                    d2 = date_start[np.where(DATE.day == date_start)[0]+1][0]-1
                    d2str = "%02d"%d2
                elif DATE.day == date_start[-1]:
                    d1str = "%02d"%DATE.day
                    d2 = monthrange(DATE.year,DATE.month)[-1]
                    d2str = "%02d"%d2
                else:
                    continue
                
                fn = '3HRLY/'+ystr+'/NARR'+file_collection_name+'_'+ystr+mstr+'_'+ \
                d1str+d2str+'.tar'
                fn_list = np.hstack((fn_list,fn))
        fn_list = np.ndarray.tolist(fn_list)
        self.tar_list = fn_list
    
    def F_download_narr(self,ret):
        """
        download narr files in filelist, a list generated from F_generate_list
        ret:
            an object generated by 
            url = 'https://rda.ucar.edu/cgi-bin/login'
            values = {'email' : login_email, 'passwd' : login_password, 'action' : 'login'}
            ret = requests.post(url,data=values)
        updated on 2019/05/29
        """
        cwd = os.getcwd()
        os.chdir(self.narr_dir)
        
        dspath = 'http://rda.ucar.edu/data/ds608.0/'
        for file in self.tar_list:
            filename=dspath+file
            file_base = os.path.basename(file)
            self.logger.info('Downloading',file_base)
            req = requests.get(filename, cookies = ret.cookies, allow_redirects=True, stream=True)
            filesize = int(req.headers['Content-length'])
            with open(file_base, 'wb') as outfile:
                chunk_size=1048576
                for chunk in req.iter_content(chunk_size=chunk_size):
                    outfile.write(chunk)
                    if chunk_size < filesize:
                        check_file_status(file_base, filesize)
                        #check_file_status(file_base, filesize)
        os.chdir(cwd)
    
    def F_untar_filelist(self,if_delete_tar=False):
        """
        untar the downloaded narr files. may delete .tar files after extraction
        if_delete_tar:
            delete tar files or not
        updated on 2019/05/29
        """
        import tarfile
        cwd = os.getcwd()
        os.chdir(self.narr_dir)
        for file in self.tar_list:
            fn = os.path.join(self.narr_dir,file.split('/')[-1])
            self.logger.info('untaring '+fn)
            tar = tarfile.open(fn)
            tar.extractall()
            tar.close()
            if if_delete_tar:
                self.logger.info('deleting '+fn)
                os.remove(fn)
        os.chdir(cwd)    
    
    def F_load_narr(self,file_collection_names=['sfc','flx','3D'],\
                    file_collection_fields=[['HPBL_221_SFC','PRES_221_SFC','TMP_221_SFC','FRICV_221_SFC','SHTFL_221_SFC','LHTFL_221_SFC'],\
                                            ['PRES_221_TRO','HGT_221_TRO','U_GRD_221_HTGL','V_GRD_221_HTGL'],\
                                            ['HGT_221_ISBL','SPF_H_221_ISBL','TMP_221_ISBL','TKE_221_ISBL','U_GRD_221_ISBL','V_GRD_221_ISBL']],\
                    rename_fields=[['PBLH','P_surf','T_surf','u_star','sensible_heat','latent_heat'],\
                                   ['P_tropopause','GPH_tropopause','U_10m30m','V_10m30m'],\
                                   ['GPH','q','Temperature','TKE','U','V']]):
        """
        load narr files to memory and subset
        file_collection_names:
            a list of narr collection names, such as ['sfc','3D']
        file_collection_fields:
            variable names in grib files, try F_show_grib_variables(fn) to see all 
            possible variable names used in a grib file
        rename_fields:
            what you want to call those variables, have to be the same size as file_collection_fields
        """
        self.file_collection_names = file_collection_names
        self.file_collection_fields = file_collection_fields
        self.rename_fields = rename_fields
        narr_data = {}
        # allocate memory for geos_data        
        for i in range(len(file_collection_names)):
            file_collection_name = file_collection_names[i]
            file_collection_field = rename_fields[i]
            
            for ifield in range(len(file_collection_field)):
                fn = file_collection_field[ifield]
                fidx = file_collection_fields[i][ifield]
                if file_collection_name == '3D':
                    narr_data[fn] = np.zeros((self.ny,self.nx,self.nlayer,self.nstep))
                # some wind in flx are stacked, pretty annoying
                elif file_collection_name == 'flx' and (fidx == 'U_GRD_221_HTGL' or fidx == 'V_GRD_221_HTGL'):
                    narr_data[fn] = np.zeros((self.ny,self.nx,2,self.nstep))
                else:
                    narr_data[fn] = np.zeros((self.ny,self.nx,self.nstep))
        
        self.matlab_datenum = np.zeros((self.nstep),dtype=np.float64)
        
        for istep in range(self.nstep):
            file_datetime = self.narr_start_datetime+datetime.timedelta(hours=self.step_hour*istep)
            file_datenum = (file_datetime.toordinal()\
                                +file_datetime.hour/24.\
                                +file_datetime.minute/1440.\
                                +file_datetime.second/86400.+366.)
            self.matlab_datenum[istep] = file_datenum
            
            for i in range(len(file_collection_names)):
                if file_collection_names[i] == '3D':
                    file_path = os.path.join(self.narr_dir,'merged_AWIP32'+\
                                             file_datetime.strftime('.%Y%m%d%H.')+\
                                             file_collection_names[i])
                else:
                    file_path = os.path.join(self.narr_dir,'merged_AWIP32'+\
                                             file_datetime.strftime('.%Y%m%d%H.RS.')+\
                                             file_collection_names[i])
                self.logger.info('loading '+file_path)
                file_collection_field= file_collection_fields[i]
                rename_field = rename_fields[i]
                outp = F_open_narr_grib(file_path,file_collection_field,rename_field,if_save_coordinates=False)
                for ifield in range(len(file_collection_field)):
                    #self.logger.debug(rename_field[ifield]+' grib field: '+file_collection_field[ifield])
                    #self.logger.debug(narr_data[rename_field[ifield]].shape)
                    narr_data[rename_field[ifield]][...,istep] = F_subset_narr_variable(\
                             outp[rename_field[ifield]],\
                             self.startx_index_1based,self.endx_index_1based,\
                             self.starty_index_1based,self.endy_index_1based)
        
        self.narr_data = narr_data
    
    def F_save_narr_data2mat(self,if_delete_grib=False):
        """
        save narr_data loaded by F_load_narr to mat file
        need some special treatment for U/V_10m30m
        created on 2019/05/30
        """
        if self.nstep != 1:
            self.logger.error('nstep = '+'%d'%self.nstep+', this function only works for single time step (start=end)!')
            return
        file_datetime = self.narr_start_datetime
        file_dir = os.path.join(self.narr_dir,file_datetime.strftime('Y%Y'),\
                                   file_datetime.strftime('M%m'),\
                                   file_datetime.strftime('D%d'))
        if not os.path.exists(file_dir):
            self.logger.warning('directory '+file_dir+' does not exist! creating one...')
            os.makedirs(file_dir)
        mat_fn = os.path.join(file_dir,'subset_'+file_datetime.strftime('%Y%m%d_%H%M')+'.mat')
        save_dict = self.narr_data
        save_dict = {k:(np.squeeze(v)).astype(np.float32) for (k,v) in save_dict.items()}
        save_dict['x'] = self.x.flatten()
        save_dict['y'] = self.y.flatten()
        save_dict['x0'] = self.x0.flatten()
        save_dict['y0'] = self.y0.flatten()
        save_dict['startx_index_1based'] = self.startx_index_1based
        save_dict['endx_index_1based'] = self.endx_index_1based
        save_dict['starty_index_1based'] = self.starty_index_1based
        save_dict['endy_index_1based'] = self.endy_index_1based
        if 'U_10m30m' in save_dict.keys():
            save_dict['U_10m'] = save_dict['U_10m30m'][...,0]
            save_dict['U_30m'] = save_dict['U_10m30m'][...,1]
            del save_dict['U_10m30m']
            #self.logger.debug('U_10m30m separated and deleted')
        if 'V_10m30m' in save_dict.keys():
            save_dict['V_10m'] = save_dict['V_10m30m'][...,0]
            save_dict['V_30m'] = save_dict['V_10m30m'][...,1]
            del save_dict['V_10m30m']
            #self.logger.debug('V_10m30m separated and deleted')
        sio.savemat(mat_fn,save_dict)
        if not if_delete_grib:
            return
        for i in range(len(self.file_collection_names)):
            if self.file_collection_names[i] == '3D':
                self.file_path = os.path.join(self.narr_dir,'merged_AWIP32'+\
                                         file_datetime.strftime('.%Y%m%d%H.')+\
                                         self.file_collection_names[i])
            else:
                file_path = os.path.join(self.narr_dir,'merged_AWIP32'+\
                                         file_datetime.strftime('.%Y%m%d%H.RS.')+\
                                         self.file_collection_names[i])
            self.logger.info('deleting '+file_path)
            os.remove(file_path)
    