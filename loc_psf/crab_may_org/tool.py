import os
import glob

def get_rawdata(data_dir, instrument="HE"):
    if instrument == "HE":
        try:filename     = sorted(glob.glob(data_dir + '/HE/*HE-Evt_FFFFFF_V[1-9]*'))[-1]
        except:print("\nERROR: Event file(Evt) not exist...skip this observation\n");
        try:orbitname    = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
        except:print("\nERROR: Orbit file(Orbit) not exist...skip this observation\n");
        try:attname      = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
        except:print("\nERROR: Attitude file(Att) not exist...skip this observation\n");
        try:hvfilename   = sorted(glob.glob(data_dir + '/HE/HXMT*HV_FFFFFF*V[1-9]*'))[-1]
        except:print("\nERROR: High Voltage file(HV) not exist...skip this observation\n");
        try:pmfilename   = sorted(glob.glob(data_dir + '/HE/HXMT*PM_FFFFFF*V[1-9]*'))[-1]
        except:print("\nERROR: Particle Monitor file(PM) not exist...skip this observation\n");
        try:deadfilename = sorted(glob.glob(data_dir + '/HE/HXMT*DTime*V[1-9]*'))[-1]
        except:print("\nERROR: Dead time file(DTime) not exist...skip this observation\n");
        try:tempfilename = sorted(glob.glob(data_dir + '/HE/HXMT*TH*V[1-9]*'))[-1]
        except:print("\nERROR: Temperature file(TH) not exist...skip this observation\n");
        try:ehkfilename  = sorted(glob.glob(data_dir + '/AUX/*_EHK_*V[1-9]*'))[-1]
        except:print("\nERROR: House Keeping file(EHK) not exist...skip this observation\n");
        rawdata_name = ["EVT", "Orbit", "ATT", "HV", "PM", "DTime", "TH", "EHK"]
        rawdata_content = [filename, orbitname, attname, hvfilename, pmfilename, deadfilename, tempfilename, ehkfilename]
        rawdata_dict = dict(zip(rawdata_name, rawdata_content))
        return rawdata_dict
    elif instrument == "ME":
        try:filename     = sorted(glob.glob(data_dir + '/ME/*ME-Evt_FFFFFF_V[1-9]*'))[-1]
        except:print("\nERROR: Event file(Evt) not exist...skip this observation\n");
        try:orbitname    = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
        except:print("\nERROR: Orbit file(Orbit) not exist...skip this observation\n");
        try:attname      = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
        except:print("\nERROR: Attitude file(Att) not exist...skip this observation\n");
        try:tempfilename = sorted(glob.glob(data_dir + '/ME/HXMT*TH*V[1-9]*'))[-1]
        except:print("\nERROR: Temperature file(TH) not exist...skip this observation\n");
        try:instatname   = sorted(glob.glob(data_dir + '/ME/HXMT*InsStat*V[1-9]*'))[-1]
        except:print("\nERROR: Instrument Status file(InsStat) not exist...skip this observation\n");
        try:ehkfilename  = sorted(glob.glob(data_dir + '/AUX/*_EHK_*V[1-9]*'))[-1]
        except:print("\nERROR: House Keeping file(EHK) not exist...skip this observation\n");
        rawdata_name = ["EVT", "Orbit", "ATT", "TH", "EHK", "InsStat"]
        rawdata_content = [filename, orbitname, attname, tempfilename, ehkfilename, instatname]
        rawdata_dict = dict(zip(rawdata_name, rawdata_content))
        return rawdata_dict
    elif instrument == "LE":
        try:filename     = sorted(glob.glob(data_dir + '/LE/*LE-Evt_FFFFFF_V[1-9]*'))[-1]
        except:print("\nERROR: Event file(Evt) not exist...skip this observation\n");
        try:orbitname    = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
        except:print("\nERROR: Orbit file(Orbit) not exist...skip this observation\n");
        try:attname      = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
        except:print("\nERROR: Attitude file(Att) not exist...skip this observation\n");
        try:tempfilename = sorted(glob.glob(data_dir + '/LE/HXMT*TH*V[1-9]*'))[-1]
        except:print("\nERROR: Temperature file(TH) not exist...skip this observation\n");
        try:instatname   = sorted(glob.glob(data_dir + '/LE/HXMT*InsStat*V[1-9]*'))[-1]
        except:print("\nERROR: Instrument Status file(InsStat) not exist...skip this observation\n");
        try:ehkfilename  = sorted(glob.glob(data_dir + '/AUX/*_EHK_*V[1-9]*'))[-1]
        except:print("\nERROR: House Keeping file(EHK) not exist...skip this observation\n");
        rawdata_name = ["EVT", "Orbit", "ATT", "TH", "EHK", "InsStat"]
        rawdata_content = [filename, orbitname, attname, tempfilename, ehkfilename, instatname]
        rawdata_dict = dict(zip(rawdata_name, rawdata_content))
        return rawdata_dict


