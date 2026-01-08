class Log:
    """
    Centralized logging utility for CEDAR.

    This class provides simple hierarchical logging with optional
    console output and file logging. Log indentation levels can be
    increased or decreased to reflect nested operations (e.g.,
    model setup, solver iterations).

    Logging is performed entirely through class methods and shared
    class-level state.

    Attributes
    ----------
    problem_path : str
        Root path for the current problem.

    levels : int
        Current indentation level used for hierarchical logging.

    log_to_file : bool
        If True, log messages are written to a log file.

    print_messages : bool
        If True, log messages are printed to stdout.

    is_start : bool
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
          ||||          Computational Environment for Dynamics
          ||||          Analysis of Reactor systems in space

          
"""

    problem_path = ""
    levels = 0
    log_to_file = False
    print_messages = False
    is_start = True

    @classmethod
    def create(cls, problem_path: str):
        """
        Initialize logging for a problem.

        Enables file logging and console output and creates the
        log file at ``<problem_path>/outputs/cedar.log``.

        Parameters
        ----------
        problem_path : str
            Root directory of the problem.
        """
        cls.log_to_file = True
        cls.print_messages = True
        cls.log_path = problem_path + "/outputs/cedar.log"

        cls._create_log_file()

    @classmethod
    def add_level(cls):
        """
        Increase the current logging indentation level.

        Typically used when entering a nested operation.
        """
        cls.levels += 1

    @classmethod
    def remove_level(cls):
        """
        Decrease the current logging indentation level.

        Typically used when exiting a nested operation.
        """

        if cls.levels == 0:
            return

        cls.levels -= 1

    @classmethod
    def message(cls, message: str):
        """
        Write a log message.

        The message is indented according to the current logging
        level and written to stdout and/or the log file depending
        on configuration.

        Parameters
        ----------
        message : str
            Message to log.
        """
        if cls.is_start:
            cls.is_start = False
            cls.message(cls.start_message)

        if cls.levels != 0:
            message = "       " * cls.levels + message

        if cls.print_messages:
            print(message)

        cls._add_to_log_file(message)

    @classmethod
    def line_break(cls):
        """
        Insert a blank line in the log output.
        """
        if cls.print_messages:
            print("")

        cls._add_to_log_file("")

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
        cls._add_to_log_file(message)

        raise Exception(message)

    @classmethod
    def _create_log_file(cls):
        """
        Create or overwrite the log file.
        """
        with open(cls.log_path, "w") as file:
            file.write("")

    @classmethod
    def _add_to_log_file(cls, message: str):
        """
        Append a message to the log file if file logging is enabled.

        Parameters
        ----------
        message : str
            Message to append.
        """
        if cls.log_to_file:
            with open(cls.log_path, "a") as file:
                file.write(message + "\n")