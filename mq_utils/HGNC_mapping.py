#!/usr/bin/python3


import httplib2 as http
import json
import pandas as pd
from ratelimit import limits, sleep_and_retry

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

# 10 calls per second
CALLS = 10
RATE_LIMIT = 1


@sleep_and_retry
@limits(calls=CALLS, period=RATE_LIMIT)
def get_HGNC_mapping(id, request="symbol"):
    headers = {
        'Accept': 'application/json',
    }

    output = pd.DataFrame()

    uri = 'http://rest.genenames.org'
    path = '/fetch/' + request + '/' + id

    target = urlparse(uri + path)
    method = 'GET'
    body = ''

    h = http.Http()

    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)

    if response['status'] == '200':
        status = "success"
    else:
        print('Error detected: ' + response['status'])
        status = "error"
        return None

    if status == "success":
        data = json.loads(content)

        if (len(data["response"]["docs"]) == 0):
            return None
        else:
            for entry in data["response"]["docs"]:
                hgnc_id = entry.get("hgnc_id", None)
                symbol = entry.get("symbol", None)
                prev = entry.get("prev_symbol", None)
                alias = entry.get("alias_symbol", None)

                if prev is not None:
                    prev = ";".join(prev)

                if alias is not None:
                    alias = ";".join(alias)

                output_chunk = pd.Series([request, id, hgnc_id, symbol, prev, alias])
                output = pd.concat([output, output_chunk], axis=1)

            output = output.T
            output.columns = ["Request Type", "Input", "HGNC ID", "Symbol", "Previous Symbol", "Alias Symbol"]
            return (output)
