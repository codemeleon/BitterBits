#!/usr/bin/env python
"""Find commands."""
print("Anmol")
import click
import sqlite3
from os import path, makedirs
import subprocess
# import expanduser


def helper(program, kind):
    """Return help information."""
    info = ""
    if kind == "man":
        # man -P cat command_name
        info = subprocess.Popen([kind, "-P", "cat", program],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE
                                ).communicate()[0]
    else:
        info = subprocess.Popen([program, kind],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE
                                ).communicate()[0]
    return info


# TODO:60 need to add regular expression and multisearch

@click.command()
@click.option("-upgrade", help="Upgrade the database", type=bool,
              default=False, show_default=True)
@click.option("-cmdonly", help="Show command details", type=bool,
              default=True, show_default=True)
@click.option("-search", help="Search term", type=str,
              default=None, show_default=True)
def run(upgrade, cmdonly, search):
    """Generate table for coammand and allows search based on words."""
    home = path.expanduser("~")
    if not path.isdir("%s/.whatapp" % home) or upgrade:
        if not path.isdir("%s/.whatapp" % home):
            makedirs("%s/.whatapp" % home)
        conn = sqlite3.connect('%s/.whatapp/whatapp.sqlite3' % home)
        cursor = conn.cursor()
        if not path.isdir("%s/.whatapp" % home):
            cursor.execute("create table commands (command text,\
                           details text)")
        else:
            cursor.execute("drop table if exists commands")
            cursor.execute("create table commands (command text,\
                           details text)")
        click.echo("Populating/updating data table...")
        programs_all = subprocess.check_output("compgen -c", shell=True,
                                               executable="bash")
        programs_all = programs_all.splitlines()
        programs_all = [x.decode('ascii') for x in programs_all]
        programs_alias = subprocess.check_output("compgen -a", shell=True,
                                                 executable="bash")
        programs_alias = programs_alias.splitlines()
        programs_alias = [x.decode('ascii') for x in programs_alias]
        programs = set(programs_all) - set(programs_alias)
        t_programs = []
        for program in programs:
            if len(program) < 2 and program != "R":
                continue
            if (program.startwith("_") or
                    program.startwith("{") or
                    program.startwith("}") or
                    program.startwith("[") or
                    program.startwith("]") or
                    program.startwith(".") or
                    program.startwith("!")):
                    continue
            t_programs.append(program)

        for program in t_programs:
            info = helper(program, "man")
            if len(info) < 200:
                info = helper(program, "--help")
                if len(info) < 200:
                    info = helper(program, "-help")
                    if len(info) < 200:
                        info = helper(program, "--h")
                        if len(info) < 200:
                            info = helper(program, "-h")
                            if len(info) < 200:
                                info = None
            info = info.splitlines()
            info = "\n".join(info)
            cursor.execute("INSERT INTO commands values (program, info)")
        conn.commit()
        conn.close()
        click.echo("Please enter your query again.")

    # conn = sqlite3.connect('%s/.whatapp/whatapp.sqlite3' % home)
    # cursor = conn.cursor()
    # cursor.execute("")
    # conn.close()
    # Search the table and produce the result


if __name__ == '__main__':
    run()
