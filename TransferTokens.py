import os
import argparse
import subprocess


def readIpAddresses(f):
    '''
    Description: Read host file of Ansible, fetch ip-addresses listed under [reseda]
    Input: file name of the hosts file
    Output: list of ip addresses
    '''
    start_reading = 0
    ip_addresses = list()

    # Read host file of Ansible, fetch ip-addresses listed under [reseda]
    fh = open(f)
    for line in fh:
        line = line.rstrip() # remove newline

        if line == "": # skip empty lines
            continue
        if line.startswith("#"): # skip lines that are commented out
            continue

        # start reading after encountering this keyword
        if line == "[reseda]":
            start_reading = 1
            continue

        # if there is another category, then stop reading
        if start_reading == 1 and line.startswith("["):
            break

        # Store ip addresses
        if start_reading == 1:
            if "#" in line: # remove comments after the ip-address
                line = line.split("#")[0]
                line = line.rstrip()
            ip_addresses.append(line)
    fh.close()

    return(ip_addresses)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Transfers tokens (jobs) to cloud machines')
    parser.add_argument('-d', '--directory', default='tokens/', type=str, help='Directory with tokens, files start with SAMPLE- (default: %(default)s)')
    parser.add_argument('-i', '--inventory', default='../ansible-playbooks/Reseda/hosts', type=str, help='File containing ip-addresses, also used for Ansible (default: %(default)s)')
    args = parser.parse_args()

    # Read list of ip addresses for the Ansible hosts file
    ip_addresses = readIpAddresses(args.inventory)

    # Fetch file names of the tokens that need to be transferrred to the cloud machines
    mytokens = [x for x in os.listdir(args.directory) if x.startswith("SAMPLE")]

    ## Divide the jobs/tokens over the ip addresses

    # Initialize variable to contain the hosts and ip-addresses
    ip_tokens = dict()
    for ip in ip_addresses:
        ip_tokens[ip] = list()

    # Store tokens per ip address
    nr_hosts = len(ip_addresses)
    for token_index in range(len(mytokens)):
        host_index = token_index % nr_hosts
        ip_tokens[ip_addresses[host_index]].append(args.directory + mytokens[token_index])

    # Transfer the tokens to the cloud machines
    for ip in ip_tokens:
        command = "scp " + " ".join(ip_tokens[ip]) + " bschaikvan@" + ip + ":~/git/reseda/tokens/"
        print(command)
        command = command.split(" ")
        subprocess.run(command)

