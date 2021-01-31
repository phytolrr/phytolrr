import logging
import urllib.request

ipv4_url = 'https://www.cloudflare.com/ips-v4'
ipv6_url = 'https://www.cloudflare.com/ips-v6'


def get_ips(url):
    headers = {'User-Agent': 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'}
    req = urllib.request.Request(url, headers=headers)
    resp = urllib.request.urlopen(req)
    if resp.getcode() != 200:
        logging.error("Failed to get ips from url " + url)
        return []
    ips_bytes = resp.read()
    if len(ips_bytes) == 0:
        logging.error(str.format("Failed to get ips from url {}, no content", url))
        return []
    ips_str = ips_bytes.decode('utf8')
    if len(ips_str) == 0:
        logging.error(str.format("Failed to get ips from url {}, invalid content {}", url, ips_bytes))
    ips = ips_str.splitlines()
    ips.sort()
    return ips


ips_v4 = get_ips(ipv4_url)
ips_v6 = get_ips(ipv6_url)


FILE_NAME = 'real_ip_cloudflare.conf'

with open(FILE_NAME, 'w') as f:
    f.writelines(['set_real_ip_from ' + ip + ';\n' for ip in ips_v4])
    f.writelines(['set_real_ip_from ' + ip + ';\n' for ip in ips_v6])
    f.writelines([
        '\n',
        '# use any of the following two\n',
        'real_ip_header CF-Connecting-IP;\n',
        '#real_ip_header X-Forwarded-For;\n'
    ])
