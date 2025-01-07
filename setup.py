
def get_long():
    import os
    with open(
        os.path.join(os.path.dirname(__file__), "README.rst"), encoding="utf-8"
    ) as f:
        readme = f.read()
    with open(
        os.path.join(os.path.dirname(__file__), "CHANGELOG.rst"), encoding="utf-8"
    ) as f:
        change_log = f.read()
    return "%s\n%s" % (readme, change_log)

def get_requirements():
    import os
    with open(
        os.path.join(os.path.dirname(__file__), "requirements.txt"), encoding="utf-8"
    ) as f:
        return list(f.read().splitlines())

def setup():
    import setuptools, os, glob
    setuptools.setup(
        name="idrfeatlib",
        version="0.0.0",
        description="A hopefully usable internal IDR feature analysis library.",
        long_description=get_long(),
        long_description_content_type="text/x-rst",
        author="Tian Hao Huang",
        author_email="tianh.huang@mail.utoronto.ca",
        url="https://github.com/alex-tianhuang/usable_featlib",
        packages=setuptools.find_packages("src"),
        package_dir={"": "src"},
        py_modules=[os.path.splitext(os.path.basename(i))[0] for i in glob.glob("src/*.py")],
        include_package_data=True,
        zip_safe=False,
        install_requires=get_requirements(),
    )

setup()