#!/usr/bin/expect -f

#process to parse command line arguments into OPTS array
proc parse_args {argc argv} {
    global OPTS
    foreach {key val} $argv {
        if {$key eq "-u" || $key eq "--unset"} {
            clear_creds
            exit 0
        }
    }
}

# Checks if system is macOS
proc is_macos {} {
    set result [catch {exec echo $OSTYPE} output]
    if {$result eq "darwin*"} {
        return true
    } else {
        return false
    }
}

# Fetches existing env variable
proc get_env_var {var_name} {
    set value [exec sh -c "echo \$$var_name"]
    return $value
}

# Helper for if var defined
proc has_pos_length {value} {
    if {[string length $value] == 0} {
        return false
    } else {
        return true
    }
}

# Prompts user for var if not defined
proc prompt_env_var {var_name prompt} {
    set value [get_env_var $var_name]
    if {[has_pos_length $value] == false} {
        send_user "$prompt"
        flush stdout
        gets stdin value

        # Save to env
        save_to_env $var_name $value
    }
    return $value
}

# Gets the current shell (e.g. '/usr/bin/zsh')
proc get_shell {} {
    set shell [exec sh -c "echo \$SHELL"]
    return $shell
}

# Get shell name (e.g. 'zsh')
proc get_shell_name {} {
   set shell [get_shell]   
   regexp {[^/]+$} $shell shell_name

   return $shell_name
}

# Gets the current config file (e.g. ~/.zshrc)
proc get_config_file {} {
    set shell_name [get_shell_name]
    set config_file ""

    # Determine the configuration file based on the detected shell
    set config_file "~/.${shell_name}rc"

    if {[is_macos] == true} {
        set config_file "~/.bash_profile"
    }
    return $config_file
}

# Saves an env variable to shell run commands
proc save_to_env {var_name value} {
    set config_file [get_config_file]
    
    # Remove this variable incase it's already set
    rmv_env_from_rc $var_name false
   
    # Add to config
    exec sh -c "echo 'export $var_name=\"$value\"' >> $config_file"
    send_user "Saved $var_name to $config_file\n"
}

# Reload shell config to apply changes made to run commands
proc reload_config {} {
    send_user "Reloading config...\n"
    set config_file [get_config_file]

    send_user "Manually reload your config via '. ${config_file}'\n\n"
}

# Removes and env var from the run command profile (e.g. ~/.bashrc)
proc rmv_env_from_rc {var_name disp_log} {
    set config_file [get_config_file]

    if {$disp_log eq true} {
        send_user "    Removing $var_name from $config_file\n"
    }
    
    set var_exists [exec sh -c "grep -q '^export $var_name=' $config_file; echo \$?"]
    if {$var_exists == 0} {
        exec sh -c "sed -i '/^export $var_name=/d' $config_file"
        if {$disp_log eq true} {
            send_user "Cleared ${var_name} from ${config_file}\n"
        }
    } 
}

# Clears username and password from env
proc clear_creds {} {
    send_user "Clearing credentials from run commands...\n"
    rmv_env_from_rc "GETAFIX_USR" true
    rmv_env_from_rc "GETAFIX_PWD" true
    
    send_user "Note you will still need to run the following commands to remove the credentials from your env vars\n\n"
    
    send_user "   - unset GETAFIX_USR GETAFIX_PWD\n\n"

    send_user "You can verify this has worked by running
   - printenv | grep GETAFIX_

...and confirming no results appear.\n\n"
}

# Masks a password
proc mask_pwd {pass} {
    set pwd_length [string length $pass]
    set masked_pwd [string repeat "*" $pwd_length]
    return $masked_pwd
}

# Entry

# Parse command line args
parse_args $argc $argv

set timeout -1

# Get username and password from env
set has_usr [has_pos_length [get_env_var "GETAFIX_USR"]]
set has_pwd [has_pos_length [get_env_var "GETAFIX_PWD"]]

set user [prompt_env_var "GETAFIX_USR" "Username not set. Please enter your username: "]
set password [prompt_env_var "GETAFIX_PWD" "Password not set. Please enter your password: "]

if {$has_usr == false || $has_pwd == false} {
    # Need to reload shell
    reload_config
    exit 0
}

set host "getafix.smp.uq.edu.au"

send_user "SSH Config:\n"
send_user "    - Username: $user\n"
send_user "    - Password: [mask_pwd $password]\n"

send_user "\n"
send_user "Connecting to $user@$host..."

# Enable logging of commands to and from remote
log_user 1;

# Spawn ssh session
spawn ssh -oStrictHostKeyChecking=no -oCheckHostIP=no $user@$host

expect {
	-re "(P|p)assword: *$" {
		send -- "$password\r"
        interact;
	}
	timeout {
		send_user "Error: Connection timed out.\n"
		exit 1
	}
	eof  {
		send_user "Error Unexpected end of file.\n"
		exit 1
	}
}



