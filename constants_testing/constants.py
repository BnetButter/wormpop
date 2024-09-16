import json

# Global variable to hold the path to the config file
CONFIG_PATH = None
_config = None  # Cached config

def set_config_path(path):
    """Sets the global config path."""
    global CONFIG_PATH
    CONFIG_PATH = path

def load_constants():
    """Loads the constants from the set config path."""
    global _config

    if _config is None:  # Load the config only once
        if CONFIG_PATH is None:
            raise ValueError("Config path is not set. Please call 'set_config_path()' first.")
        
        with open(CONFIG_PATH, 'r') as file:
            _config = json.load(file)
    
    return _config