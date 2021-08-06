import chardet
from os import walk
from os.path import join


def main_subdirs():
    subdirs = set()
    for subdir, dirs, files in walk("maildir"):
        if "/" in subdir:
            subdirs.add(subdir.split("/")[1])
    return subdirs

def get_relevant_enron_email_addresses():
    emails_for_employees = []  # (employee, email_address_candidates)
    subdirs = sorted(list(main_subdirs()))
    print("Got %d subdirs" % len(subdirs))
    for subdir in subdirs:  # Unique employees
        last_name_s = subdir.split("-")[:-1]
        found_files_in_sent = False
        candidate_email_counts = {}
        for the_dir, _, files in walk("maildir/" + subdir):
            if ("sent" not in the_dir) and "deleted_items" not in the_dir:
                # Someone, I forget who, has no sent directories, but they do
                #   have sent emails in the deleted_items folder.
                continue
            for filename in files:
                found_files_in_sent = True
                full_file_path = join(the_dir, filename)
                # First, find the filetype
                with open(full_file_path, "rb") as f:
                    filetype = chardet.detect(f.read())

                f = open(full_file_path, "r", encoding=filetype["encoding"])
                lines = f.readlines()
                f.close()

                sender = lines[2].split(" ")[-1]
                sender = strip_to_email_address(sender)
                if "@enron.com" not in sender:
                    continue

                has_a_last_name = False
                for last_name in last_name_s:
                    if last_name == "ybarbo":
                        last_name = "y'barbo"
                    last_name_no_apost = last_name.replace("'", "")
                    if last_name_no_apost in sender:
                        sender = sender.replace(last_name_no_apost, last_name)
                        has_a_last_name = True
                        break
                if not has_a_last_name:
                    continue

                if sender not in candidate_email_counts:
                    candidate_email_counts[sender] = 0
                candidate_email_counts[sender] += 1
        candidate_email_counts = \
            [(c, e) for e, c in candidate_email_counts.items()]
        candidate_email_counts.sort(reverse=True)
        candidate_email_counts = [(e, c) for (c, e) in candidate_email_counts]
        has_first_name = None
        has_first_initial = False
        has_last_name_only = False
        cec = []
        for (e, c) in candidate_email_counts:
            if str(e[1]) == ".":
                if not has_first_initial:
                    cec.append((e, c))
                    has_first_initial = True
            elif str(e[0]) == ".":
                if not has_last_name_only:
                    cec.append((e, c))
                    has_last_name_only = True
            else:
                first_name_in_e = e.split(".")[0]
                if has_first_name is None:
                    has_first_name = first_name_in_e
                    cec.append((e, c))
                elif (first_name_in_e in has_first_name) or \
                        (has_first_name in first_name_in_e):
                    # This case is for bert/albert meyers
                    cec.append((e, c))
        emails_for_employees.append((subdir, cec))
        print("")
        print(emails_for_employees[-1])
    return emails_for_employees


def strip_to_email_address(word):
    word = word.replace(",", "")
    word = word.replace("'", "")
    word = word.replace("<", "")
    word = word.replace(">", "")
    word = word.replace("\n", "")
    return word

def month_to_days(month, is_leap_year):
    if "jan" in month.lower():
        return 0
    elif "feb" in month.lower():  # 31
        return 31
    elif "mar" in month.lower():  # 28
        return 59 + int(is_leap_year)
    elif "apr" in month.lower():  # 30
        return 89 + int(is_leap_year)
    elif "may" in month.lower():  # 31
        return 120 + int(is_leap_year)
    elif "jun" in month.lower():  # 30
        return 150 + int(is_leap_year)
    elif "jul" in month.lower():  # 31
        return 181 + int(is_leap_year)
    elif "aug" in month.lower():  # 31
        return 212 + int(is_leap_year)
    elif "sep" in month.lower():  # 30
        return 242 + int(is_leap_year)
    elif "oct" in month.lower():  # 31
        return 273 + int(is_leap_year)
    elif "nov" in month.lower():  # 30
        return 303 + int(is_leap_year)
    elif "dec" in month.lower():  # 31
        return 334 + int(is_leap_year)
    else:
        raise ValueError("Error! Unknown month: %s" % month)

# An imperfect calculation that works fine for this period.
#   It is imperfect because it makes the following assumptions:
#       1. Leap years occur every 4 years.
#       2. Year 0 was a year, but the first leap year was year 4, not year 0.
#
# Only takes PST and PDT times. (-0800 and -0700)
def date_to_seconds(date):
    year = int(date[2])
    is_leap_year = (year % 4) == 0
    month_days = month_to_days(date[1], is_leap_year)
    day = int(date[0]) - 1
    time = date[3].split(":")
    hour = int(time[0])
    if date[4] == "-0800":
        hour_offset = 0
    elif date[4] == "-0700":
        hour_offset = -1
    else:
        raise ValueError("Error! Unknown timezone offset %s" % date[4])
    minute = int(time[1])
    second = int(time[2])
    return ((((int(year * 365.25)) + month_days + day) * 24 \
                + hour + hour_offset) * 60 + minute) * 60 + second

def load_emails_to_ids():
    f = open("email_address_labels.txt", "r")
    lines = f.readlines()
    f.close()
    d = {}
    for line in lines:
        line = line.split(" ")
        i = int(line[1].strip())
        email = strip_to_email_address(line[0])
        d[email] = i
    return d

# Collects two lists of emails.
# The first list ("exclusive") has only emails where the sender and EVERY
#   recipient are in the main emails collection.
# The second list ("all") has emails where the sender and SOME of the
#   recipients are in the main emails collection. Only recipients in the main
#   emails collection are recorded.
def collect_emails():
    emails_to_ids = load_emails_to_ids()

    earliest_seconds_all = None
    first_timestamp_all = None
    earliest_seconds_exclusive = None
    first_timestamp_exclusive = None

    emails_all = set()
    emails_exclusive = set()
    message_ids = set()

    visited_names = set()

    for subdir, dirs, files in walk("maildir"):
        for filename in files:
            name = subdir.split("/")[1]
            if name not in visited_names:
                print(("Before looking at %s we have " % name) + \
                 ("%d excl, %d all" % (len(emails_exclusive), len(emails_all))))
                visited_names.add(name)
            full_file_path = join(subdir, filename)
            # if len(subdir.split("/")) < 3:
            #     continue
            # if "sent" not in subdir.split("/")[2]:
            #     continue

            # First, find the filetype
            with open(full_file_path, "rb") as f:
                filetype = chardet.detect(f.read())

            f = open(full_file_path, "r", encoding=filetype["encoding"])
            lines = f.readlines()
            f.close()

            sender = strip_to_email_address(lines[2].split(" ")[-1])
            if sender not in emails_to_ids:
                continue

            sender = emails_to_ids[sender]

            message_id = lines[0].strip()

            timestamp = lines[1].split()[2:-1]
            timestamp_in_seconds = date_to_seconds(timestamp)

            recipients = set()
            num_recipients = 0
            num_known_recipients = 0
            started_to = False
            started_cc = False
            started_bcc = False

            line_has_emails = False  # whether the current line contains emails
            for i in range(3, len(lines)):
                if len(lines[i]) >= 3 and lines[i][:3].lower() == "to:" and \
                        not started_to:
                    started_to = True
                    line_has_emails = True
                elif len(lines[i]) >= 3 and lines[i][:3].lower() == "cc:" and \
                        not started_cc:
                    started_cc = True
                    line_has_emails = True
                elif len(lines[i]) >= 4 and lines[i][:4].lower() == "bcc:" and \
                        not started_bcc:
                    started_bcc = True
                    line_has_emails = True
                elif len(lines[i]) > 1 and str(lines[i][0]) == " ":
                    line_has_emails = line_has_emails  # continue as before
                else:
                    line_has_emails = False

                if not line_has_emails:
                    continue

                for word in lines[i].split(" "):
                    word = strip_to_email_address(word)
                    if "@" in word:
                        num_recipients += 1
                    if word in emails_to_ids:
                        num_known_recipients += 1
                        recipients.add(emails_to_ids[word])

            if num_known_recipients == 0:
                continue

            if message_id in message_ids:
                print("Skipping repeat email.")
                continue
            message_ids.add(message_id)

            the_email = (timestamp_in_seconds, sender, \
                         tuple(sorted(list(recipients))))
            emails_all.add(the_email)
            if earliest_seconds_all is None or \
                    timestamp_in_seconds < earliest_seconds_all:
                earliest_seconds_all = timestamp_in_seconds
                first_timestamp_all = timestamp
            if num_recipients == num_known_recipients:
                emails_exclusive.add(the_email)
                if earliest_seconds_exclusive is None or \
                        timestamp_in_seconds < earliest_seconds_exclusive:
                    earliest_seconds_exclusive = timestamp_in_seconds
                    first_timestamp_exclusive = timestamp

    emails_all = sorted(list(emails_all))
    emails_exclusive = sorted(list(emails_exclusive))
    return (first_timestamp_exclusive, emails_exclusive, first_timestamp_all, emails_all)

if __name__ == "__main__":

    # First, get a mapping from email addresses to ID's.
    e = get_relevant_enron_email_addresses()
    f = open("email_address_labels.txt", "w")
    already_written_addresses = set()
    label = 0
    for i in range(0, len(e)):
        wrote_an_address = False
        for j in range(0, len(e[i][1])):
            (address, count) = e[i][1][j]
            if address not in already_written_addresses:
                already_written_addresses.add(address)
                f.write(address + " " + str(label))
                if j < len(e[i][1]) - 1 or i < len(e) - 1:
                    f.write("\n")
                wrote_an_address = True
            else:
                print("Duplicate! %s" % address)
        if wrote_an_address:
            label += 1
    f.close()

    # fet -- first "exclusive" timestamp
    # ee -- "exclusive" emails
    # fat -- first "all" timestamp
    # ae -- "all" emails
    (fet, ee, fat, ae) = collect_emails()
    print("First Timestamp for Exclusive Emails == %s\n" % fet)
    f = open("exclusive_emails.txt", "w")
    for i in range(0, len(ee)):
        (timestamp, sender, recipients) = ee[i]
        f.write(str(sender) + " ")
        for r in recipients:
            f.write(str(r) + " ")
        f.write(str(timestamp))
        if i < len(ee) - 1:
            f.write("\n")
    f.close()

    print("First Timestamp for Less Exclusive Emails == %s\n" % fat)
    f = open("less_exclusive_emails.txt", "w")
    for i in range(0, len(ae)):
        (timestamp, sender, recipients) = ae[i]
        f.write(str(sender) + " ")
        for r in recipients:
            f.write(str(r) + " ")
        f.write(str(timestamp))
        if i < len(ae) - 1:
            f.write("\n")
    f.close()

    """
    # Check for repeated emails with slightly different timestamps.
    f = open("less_exclusive_emails.txt", "r")
    lines = f.readlines()
    f.close()
    lines = lines[1:]
    lines = [[int(i) for i in l.strip().split(" ")] for l in lines]
    lines = [[tuple(l[:-1]), l[-1]] for l in lines]
    seconds_window = 10
    total_repeats = 0
    total_shared_senders = 0
    for i in range(0, len(lines)):
        time = lines[i][1]
        email = lines[i][0]
        for j in range(0, i):
            idx = i - (j + 1)
            t_older = lines[idx][1]
            if abs(t_older - time) > seconds_window:
                break
            email_older = lines[idx][0]
            if email == email_older:
                print("Found an email within %d seconds! :(" % seconds_window)
                print("      " + str((email_older, t_older)))
                print("  vs. " + str((email, time)))
                total_repeats += 1
                total_shared_senders += 1
                break
            elif email[0] == email_older[0]:
                # print("Found a repeat sender within %d seconds! :(" % seconds_window)
                # print("      " + str((email_older, t_older)))
                # print("  vs. " + str((email, time)))
                total_shared_senders += 1
    print("\n\nThere were a total of %d repeat emails." % total_repeats)
    print("Also, there were a total of %d repeat senders." % total_shared_senders)
    """
