import os.path
import csv
import requests
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.common.keys import Keys


def scrape_data_igem_teams(year):
    url = 'https://igem.org/Team_List?year=' + year
    source = requests.get(url).content
    soup = BeautifulSoup(source, "lxml")
    teams_table = soup.find('tbody')
    team_rows = teams_table.find_all('tr')

    # Check if file already exists
    filename = 'sponsor_scraping\\data\\team_info.csv'
    file_exists = os.path.isfile(filename)

    with open(filename, 'a', newline='') as csvfile:
        fieldnames = ['Year', 'Team', 'Wiki',
                      'Region', 'Country', 'Track', 'Kind']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if not file_exists:
            writer.writeheader()

        for team_row in team_rows:

            team_columns = team_row.find_all('td')
            team_info = []

            for team_column in team_columns:
                if team_column.find('a') != None:
                    link = team_column.find('a')
                    # print(link)
                    try:
                        if link['class'][0] == 'team_name':
                            team_info.append(team_column.text)
                    except KeyError:
                        team_info.append(link['href'])
                else:
                    team_info.append(team_column.text)

            # Pop the empty string that was appended as first
            team_info.pop(0)

            # Only add the team if the application was accepted
            if team_info[6] == 'Accepted':
                writer.writerow({'Year': year, 'Team': team_info[0], 'Wiki': team_info[1], 'Region': team_info[2],
                                'Country': team_info[3], 'Track': team_info[4], 'Kind': team_info[5]})


def scrape_data_igem_teams_all_years():
    years = ['2020', '2019', '2018', '2017', '2016', '2015', '2014', '2013', '2012', '2011', '2010',
             '2009', '2008']

    for year in years:
        scrape_data_igem_teams(year)


def run_searches():
    PATH = 'C:\\Users\\boydc\\Projects\\university\\TUDelft_iGEM_2021\\sponsor_scraping\\chromedriver.exe'
    driver = webdriver.Chrome(PATH)

    filename = 'sponsor_scraping\\data\\team_info.csv'
    file_exists = os.path.isfile(filename)
    flag = False
    if file_exists:
        with open(filename, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            i = 0
            for row in reader:
                if flag:
                    break
                if row['Region'] == 'Europe':
                    driver.get(row['Wiki'])
                    # Open a new tab
                    driver.execute_script(
                        '''window.open("http://google.com","_blank");''')
                    # Switch back to the first tab
                    driver.switch_to.window(driver.window_handles[-1])
                    i += 1
                    if i % 10 == 0:
                        result = input('Go further? [Y/n]: ')
                        if result == 'Y':
                            num_of_tabs = 10
                            for x in range(1, num_of_tabs):
                                driver.find_element_by_tag_name(
                                    'body').send_keys(Keys.COMMAND + 'W')
                        elif result == 'n':
                            flag = True

    driver.close()


if __name__ == '__main__':
    scrape_data_igem_teams_all_years()
    # run_searches()
