import requests_html, getpass, datetime
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from bs4 import BeautifulSoup as BSoup


# function to log in to the BOLD databases, is needed to run the identification engine with >1 sequences
def bold_login(username="", password=""):
    # give user output
    print("{}: Trying to log in.".format(datetime.datetime.now().strftime("%H:%M:%S")))
    # ask for the password and username in a safe way only if it is not prvided via the input
    if not username:
        username = input("BOLD username: ")
        password = getpass.getpass("BOLD password: ")

    # start a new html session
    session = requests_html.HTMLSession()
    # update the header of the session
    session.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.0.4758.82 Safari/537.36"
        }
    )
    # define a retry strategy in case of a bad response from the BOLD server
    retry_strategy = Retry(
        total=15,
        status_forcelist=[400, 401, 403, 404, 413, 429, 500, 502, 503, 504],
        backoff_factor=1,
    )
    # create the adapter
    adapter = HTTPAdapter(max_retries=retry_strategy)
    # mount the adapter to the session
    session.mount("http://", adapter)

    # perform the post request to log in with this data
    data = {
        "name": username,
        "password": password,
        "destination": "MAS_Management_UserConsole",
        "loginType": "",
    }

    # send the post request to boldsystems.org
    session.post("http://boldsystems.org/index.php/Login", data=data)

    # test if the login was successfull
    bold_url = session.get("http://boldsystems.org")

    # parse the returned html
    soup = BSoup(bold_url.text, "html.parser")
    # look for the navigation bar
    navigation_bar = soup.find(class_="site-navigation nav navbar-nav")
    # find the text on the navigation bar
    text_fields = [field.text for field in navigation_bar.find_all("a")]
    if text_fields[5] != "Log out":
        print(
            "{}: Unable to log in.\nPlease check username and password.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
    else:
        print(
            "{}: Login successful.".format(datetime.datetime.now().strftime("%H:%M:%S"))
        )
        # return the session if the login was successful to handle the requests
        return session, username, password


if __name__ == "__main__":
    bold_login()
