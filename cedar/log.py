class Log:
    """
    Centralized logging utility for Cedar.

    Attributes
    ----------
    levels : int
        Current indentation level used for hierarchical logging.

    is_active : bool
        If True, log messages are printed to console.

    is_first_message : bool
        Indicates whether the log has just started and the start
        message has not yet been written.
    """

    start_message = """ 
           /\\            
          /**\\            _____ ______ _____          _____
         /****\\          / ____|  ____|  __ \\   /\\   |  __ \\ 
        /******\\        | |    | |__  | |  | | /  \\  | |__) |
       /********\\       | |    |  __| | |  | |/ /\\ \\ |  _  / 
      /**********\\      | |____| |____| |__| / ____ \\| | \\ \\ 
     /************\\      \\_____|______|_____/_/    \\_\\_|  \\_\\
    /**************\\   
   /****************\\  
          ||||              Computational Environment for 
          ||||              Dynamics of Advanced Reactors

"""

    levels = 0
    is_active = True
    is_first_message = True

    @classmethod
    def add_level(cls):
        """
        Increase the current logging indentation level.

        Typically used when entering a nested operation.
        """
        if cls.is_active:
            cls.levels += 1

    @classmethod
    def remove_level(cls):
        """
        Decrease the current logging indentation level.

        Typically used when exiting a nested operation.
        """
        if cls.is_active:
            cls.levels -= 1

        if cls.levels < 0:
            cls.levels = 0

        

    @classmethod
    def message(cls, message: str, end = "\n"):
        """
        Write a log message.

        Parameters
        ----------
        message : str
            Message to log.
        """
        if cls.is_first_message:
            cls.is_first_message = False
            cls.message(cls.start_message)

        if cls.levels != 0:
            message = "       " * cls.levels + message

        if cls.is_active:
            print(message, end = end, flush = (not end == "\n"))

    @classmethod
    def line_break(cls):
        """
        Insert a blank line in the log output.
        """
        if cls.is_active:
            print("")

    @classmethod
    def error(cls, message: str):
        """
        Log an error message and raise an exception.

        The message is always printed and written to the log file,
        and execution is immediately halted.

        Parameters
        ----------
        message : str
            Error message.

        Raises
        ------
        Exception
            Always raised after logging the error.
        """
        message = "ERROR :: " + message

        print(message)

        raise Exception(message)