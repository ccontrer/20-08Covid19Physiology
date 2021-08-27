
# Import libraries
import numpy as np
import pandas as pd
import platform
from urllib.request import Request, urlopen
import plotly.express as px
import dash
import dash_table as dt
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output



external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)


# Import data
req = Request('https://www.alberta.ca/data/stats/covid-19-alberta-statistics-data.csv')
req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
content = urlopen(req)
df = pd.read_csv(content)
# 'Date reported' into 'datetime' form
df['Date reported'] = pd.to_datetime(df['Date reported'], format='%Y-%m-%d')
# Fill nan values with string 'NA'
# Confirmed cases
df = df[df['Case type']=='Confirmed']

min_date = min(df['Date reported'])
max_date = max(df['Date reported'])

last_update = df['Date reported'].max()#.strftime('%B %d, %Y')
ageix = ['Under 1 year', '1-4 years', '5-9 years', '10-19 years',
     '20-29 years', '30-39 years', '40-49 years', '50-59 years',
     '60-69 years', '70-79 years', '80+ years', 'Unknown']

def plot_cases_per_date(df, ylabel):
    temp = df['Date reported'].value_counts().sort_index().fillna(0).rolling(7).mean()
    fig = px.area(temp,
                  labels = {'value':ylabel, 'index': 'Date reported'},
                  color_discrete_sequence=px.colors.qualitative.D3)
    fig.add_hline(y=500)
    fig.update_layout(showlegend=False)
    return fig

def plot_group_cases_per_date(df, group, ylabel):
    temp = df.groupby(group)['Date reported'].value_counts().unstack(0).sort_index().fillna(0).rolling(7).mean()
    if group == 'Age group':
        temp = temp.reindex(columns=ageix)
        my_color = px.colors.sequential.Viridis
    else:
        my_color = px.colors.qualitative.D3
    fig = px.area(temp,
                  labels = {'value':ylabel, 'index': 'Date reported'},
                  color_discrete_sequence=my_color)
    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01))
    return fig

# Get country names for dropdown menu (sort alphabetically)
groups = ['Alberta Health Services Zone', 'Gender', 'Age group']

app.layout = html.Div([
    # Make a title
    html.H1(
        children='Alberta COVID-19 cases dashboard',
        style={'textAlign': 'center'}
    ),
    html.H4(["Cases per day (public data available at the ", 
           html.A("Alberta website", href="https://www.alberta.ca/stats/covid-19-alberta-statistics.htm#data-export)"),
           ")."
    ]),
    html.Div(
        "Last updated: {}".format(df['Date reported'].max().strftime("%B %d, %Y"))
    ),
    # Add the plot
    html.Div([
        dcc.Graph(
            id='cases_day_graph',
            figure = plot_group_cases_per_date(df, 
                                'Gender', 
                                'Cases per day')
        )
    ]),
    # Add the plot
    html.Div([
        dcc.Graph(
            id='deaths_day_graph',
            figure = plot_group_cases_per_date(df[df['Case status']=='Died'],
                                'Gender',
                                'Deaths per day')
        )
    ]),
    html.Div([
        # Title of dropdown menu
        html.Label('Group'),
        # Format of dropdown menu
        dcc.Dropdown(id='dropdown_group',
                     options=[{'label': x, 'value': x} for x in groups],
                     value='Gender', # initial value
                     optionHeight=20,
                     searchable=False,
                     clearable=False),
    ]),
    html.Center([
        html.P("Python version: "+platform.python_version()),
        html.P("Dash version: "+dash.__version__),
    ]),
])

@app.callback(
    Output(component_id='cases_day_graph', component_property='figure'),
    [Input(component_id='dropdown_group', component_property='value')]
)
def update_figure(group):
    fig = plot_group_cases_per_date(df,
                                group,
                                'Cases per day')
    return fig

@app.callback(
    Output(component_id='deaths_day_graph', component_property='figure'),
    [Input(component_id='dropdown_group', component_property='value')]
)
def update_figure(group):
    fig = plot_group_cases_per_date(df[df['Case status']=='Died'],
                                group,
                                'Deaths per day')
    return fig

# Run app 
if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False)






