import ete3
import click


@click.command()
@click.option("-tree1", help="First input tree", type=str,
              default=None, show_default=True)
@click.option("-tree2", help="First input tree", type=str,
              default=None, show_default=True)
def run(tree1, tree2):
    assert tree1, "First tree not given. Exiting.."
    assert tree2, "Second tree not given. Exiting.."
    # Find split distance between trees
    pass

if __name__ == '__main__':
    run()
