def grab_int(linestring):
    linelist = linestring.split();
    return int(linelist[0]);

def grab_float(linestring):
    linelist = linestring.split();
    return float(linelist[0]);

def grab_string(linestring):
    linelist = linestring.split();
    return linelist[0];

def parse_ini_parameters(output_directory, 
                         parameters_file, 
                         ini_params={} ):
    """Parse a parameters.ini file.

    This function parses an input file parameters_full.ini and returns a dictionary
    containing key-value pairs for each required object in a DoGPack
    parameters.ini file.

    Input
    -----
        parameters_file - a string pointing to a valid parameters.ini file.

    Returns
    -------
        ini_params - a dictionary containing key-value pairs.  If ini_params already
        exists as a dictionary, this routine will add new keys.
    """

    import configparser;

    config = configparser.RawConfigParser();
    config.read( parameters_file );
    num_sections = len(config.sections());

    # [Params1D]
    ini_params['numframes']                = grab_int    (config['Params1D']['numframes']);
    ini_params['verbosity']                = grab_int    (config['Params1D']['verbosity']);
    ini_params['final_time']               = grab_float  (config['Params1D']['final_time']);
    ini_params['cfl']                      = grab_float  (config['Params1D']['cfl']);
    ini_params['max_number_of_time_steps'] = grab_int    (config['Params1D']['max_number_of_time_steps']);
    ini_params['Neqn']                     = grab_int    (config['Params1D']['Neqn']);
    ini_params['order']                    = grab_int    (config['Params1D']['order']);
    ini_params['Nx']                       = grab_int    (config['Params1D']['Nx']);
    ini_params['xlow']                     = grab_float  (config['Params1D']['xlow']);
    ini_params['xhigh']                    = grab_float  (config['Params1D']['xhigh']);
    ini_params['left_bc']                  = grab_string (config['Params1D']['left_bc']);
    ini_params['right_bc']                 = grab_string (config['Params1D']['right_bc']);

    return ini_params