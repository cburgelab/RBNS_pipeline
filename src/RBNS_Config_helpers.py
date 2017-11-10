#!/usr/bin/env python
import ConfigParser
import ast

def get_settings_template(
        input_config_F,
        int_settings_L,
        float_settings_L,
        int_or_float_settings_L,
        list_settings_L,
        boolean_settings_L,
        include_sections_as_first_key = False):
    """
    - A generic template that can be called by another function to return
        a settings dictionary. For example, within another function:

    def return_settings_L_D( input_config_F ):

        int_settings_L = ['num_reads_to_use'], etc. Then call:
        ...
        settings_D = get_settings_L_template(
            input_config_F,
            int_settings_L,
            float_settings_L,
            int_or_float_settings_L,
            list_settings_L,
            boolean_settings_L)
        return settings_D

    """
    Config = ConfigParser.ConfigParser()
    Config.optionxform = str
    settings_D = {}
    Config.read(input_config_F)
    for section in Config.sections():
        if (include_sections_as_first_key == True):
            settings_D[section] = {}
        options_this_section = Config.options(section)
        for option in options_this_section:
            setting = Config.get(section, option)
            if option in int_settings_L:
                setting = int(setting)

            elif option in float_settings_L:
                setting = float(setting)

            elif option in int_or_float_settings_L:
                if (float(setting) == int(float(setting))):
                    setting = int(float(setting))
                else:
                    setting = float(setting)

            elif option in list_settings_L:
                tmp_list = ast.literal_eval(setting)
                try:
                    setting = [x.strip() for x in tmp_list]
                except AttributeError:
                    setting = tmp_list
            elif option in boolean_settings_L:
                setting = Config.getboolean(section, option)
            # put it in the dictionary
            if (include_sections_as_first_key == True):
                settings_D[section][option] = setting
            else:
                settings_D[option] = setting
    return settings_D









