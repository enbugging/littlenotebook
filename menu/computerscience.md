---
layout: page
title: Computer Science
---
<ul class="posts">
    {% capture cur_year %}0{% endcapture %}
    {% for post in site.categories.computerscience %}

        {% capture year %}{{ post.date | date: '%Y' }}{% endcapture %}
        {% if year != cur_year %}
            <h3>{{ post.date | date: '%Y' }}</h3>
            {% capture cur_year %}{{year}}{% endcapture %}
        {% endif %}

        <li itemscope>
        <a href="{{ site.baseurl }}{{ post.url }}">{{ post.title }}</a>
        <p class="post-date"><span><i class="fa fa-calendar" aria-hidden="true"></i> {{ post.date | date: "%B %-d" }} - <i class="fa fa-clock-o" aria-hidden="true"></i> {% include read-time.html %}</span></p>
        </li>
    {% endfor %}
</ul>
