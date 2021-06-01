import curses


screen = curses.initscr()
curses.noecho()
screen.keypad(1)
curses.mousemask(1)
curses.cbreak()
max_y, max_x = screen.getmaxyx()
key_pressed = -2

while key_pressed != ord('q'):
    key_pressed = screen.getch()
    screen.clear()
    # if key_pressed = curses.KEY_C
    screen.addstr(0, 0, str(key_pressed))
    screen.refresh()
curses.nocbreak()
screen.keypad(0)
curses.echo()
curses.endwin()
