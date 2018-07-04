import xml.etree.ElementTree as ET
import os
import os.path
import subprocess
import tarfile

with open("SuiteSparseMatrixCollection.html") as f:
    tree = ET.parse(f)

LINKPREFIX = "https://sparse.tamu.edu/MM/"

def ensure_mtx_file(link):
    """
    Make sure the .mtx file has been downloaded. Return name of the .mtx file
    """
    tgz_name = link[len(LINKPREFIX):]
    mtx_name = tgz_name[:-len(".tar.gz")] + '.mtx' 
    if not os.path.exists(mtx_name):
        print("Downloading {}".format(tgz_name))
        os.makedirs(os.path.dirname(tgz_name), exist_ok=True)
        subprocess.run(["wget", link, '-O', tgz_name]) #    , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        with tarfile.open(tgz_name) as tar:
            for member in tar.getmembers():
                tar.extract(member, path="/tmp/")
                os.rename('/tmp/' + member.name, mtx_name)
        os.unlink(tgz_name)
    return mtx_name


def chain_decomposition(mtx_file):
    """
    Determine #cut vertex and #bridges.
    """
    cmd = "cat {} |  ~/spasm/build/bench/chains".format(mtx_file)
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
        return list(map(int, result.decode().split(";")))
    except subprocess.CalledProcessError as e:
        return -42, -42, -42

def modular_decomposition(mtx_file):
    """
    Determine #trivial modules, #nontrivial modules and average size of non-trivial modules
    """
    cmd = "cat {} |  ~/spasm/build/bench/modules".format(mtx_file)
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
        cc, trivial, nontrivial, size, largest = map(int, result.decode().split(";"))
        if nontrivial > 0:
            size /= nontrivial
        return cc, trivial, nontrivial, size, largest
    except subprocess.CalledProcessError as e:
            return -42, -42, -42, -42, -42



print("id; name; year; size; nnz; kind; CC; link; cut; bridge; CCbis; trivial; nontrivial; size; largest")
for x in tree.getroot():
    id = x.find("td[@class='column-id']").text
    name = x.find("td[@class='column-name']").text
    year = x.find("td[@class='column-date']").text
    nrows = x.find("td[@class='column-num_rows']").text
    nnz = x.find("td[@class='column-nonzeros']").text
    kind = x.find("td[@class='column-kind']").text
    links = x.findall("td/a")
    for foo in links:
        if "Matrix Market" in foo.text:
                link = foo.get('href')
    if not link.startswith(LINKPREFIX):
        raise ValueError("problème de prefixe")

    try:
        mtx_file = ensure_mtx_file(link)
        cc, cut, bridge = chain_decomposition(mtx_file)
        cc_bis, trivial, nontrivial, size, largest = modular_decomposition(mtx_file)
        #cc_bis, trivial, nontrivial, size, largest = [-1, -1, -1, -1, -1]
        print("{}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {}; {:.1f}; {}"
             .format(id, name, year, nrows, nnz, kind, link, cc, cut, bridge, cc_bis, trivial, nontrivial, size, largest))
    except ValueError:
        pass
