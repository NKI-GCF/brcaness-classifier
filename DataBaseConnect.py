import argparse,pymssql


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--server", help='enter the servern name')
parser.add_argument("-u", "--user", help='enter the username')
parser.add_argument("-p", "--pw", help='enter the password')
parser.add_argument("-d", "--db", help='enter the dbname')
parser.add_argument("-q", "--query", help='enter the query you want to execute')
args = parser.parse_args()


server = (args.server)
user = (args.user)
password = (args.pw)
database = (args.db)


conn = pymssql.connect(server, user, password,database)




cursor = conn.cursor()
cursor.execute(args.query)
conn.commit()

for row in cursor:
    print('row = %r' % (row,))

conn.close()

