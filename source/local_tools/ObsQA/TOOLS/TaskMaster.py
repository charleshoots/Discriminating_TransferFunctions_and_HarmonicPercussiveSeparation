
from modules import *
import logging,multiprocessing
import obstools
def setup_logging(log_fout):
        logging.basicConfig(filename=log_fout,level=print)
# Define the worker function
def worker(fn, args, log_file):
        # setup_logging(log_file)
        print("Worker process started with function: %s and args: %s", fn.__name__, args)
        fn(args)  # Execute the function with unpacked arguments
        print("Worker process finished.")
def main(fn,args,log_file=None):
        process = multiprocessing.Process(target=worker, args=(fn,args,log_file))
        print("Starting child process.")
        process.start()
        # Wait for the process to complete
        print("Waiting for the child process to complete.")
        process.join()  # Blocks until the child process finishes
        print('Process completed')
if __name__ == "__main__":
    # Run main program
    main()