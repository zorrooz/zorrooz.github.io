<template>
  <div class="page-section home-section">
    <header class="hero">
      <div class="hero__main">
        <span class="hero__greeting"
          ><span class="hero__greeting-mark">{{ t('greetingPrefix') }}</span
          >{{ t('greeting') }}</span
        >
        <h1 class="hero__name">{{ siteAuthor }}</h1>
        <p class="hero__bio">{{ t('developer') }}</p>
      </div>

      <div class="hero__stats">
        <div class="hero__stat">
          <span class="hero__stat-num num">{{ postCount }}</span>
          <span class="hero__stat-label">{{ t('articles') }}</span>
        </div>
        <div class="hero__stat">
          <span class="hero__stat-num num">{{ tagCount }}</span>
          <span class="hero__stat-label">{{ t('tags') }}</span>
        </div>
        <div class="hero__stat">
          <span class="hero__stat-num num">{{ totalWordsDisplay }}</span>
          <span class="hero__stat-label">{{ t('words') }}</span>
        </div>
      </div>
    </header>

    <div class="tags-row" v-if="!currentTag" v-reveal>
      <span
        v-for="tag in tagList"
        :key="tag.name"
        class="tag"
        :class="`tag--${tag.level}`"
        @click="goTag(tag.name)"
        >{{ tag.name }}</span
      >
    </div>

    <div class="posts-header" v-reveal>
      <h2 class="posts-header__title" :class="{ 'posts-header__title--tag': currentTag }">
        <template v-if="currentTag">
          <span class="posts-header__tag-title"># {{ currentTag }}</span>
        </template>
        <template v-else>{{ t('recentPosts') }}</template>
      </h2>
      <div class="posts-header__actions">
        <button
          v-if="fromPath"
          class="chip-close"
          @click="goBackFromTag"
          :aria-label="t('backToArticle')"
        >
          <i class="fas fa-arrow-left"></i>{{ t('backToArticle') }}
        </button>
        <button v-if="currentTag" class="chip-close" @click="clearTag" :aria-label="t('close')">
          <i class="fas fa-times"></i>{{ t('clearFilter') }}
        </button>
      </div>
    </div>

    <PostList :docs="filteredDocs" :perPage="5" />
  </div>
</template>

<script setup lang="ts">
defineOptions({ name: 'HomeView' })
import { computed, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useHead } from '@unhead/vue'
import { useRoute, useRouter } from 'vue-router'
import PostList from '@/components/layout/PostList.vue'
import { useLocalizedContent } from '@/composables/useLocalizedContent'
import { loadPosts, loadTags } from '@/utils/contentLoader'
import { goToTag } from '@/utils/tagQuery'
import { scrollToTop } from '@/utils/scroll'
import { SITE, toSupportedLocale } from '@/config/site'

useHead({ title: 'zorrooz’s blog - Home' })

const { t, locale } = useI18n()
const route = useRoute()
const router = useRouter()

const { data: postData } = useLocalizedContent(() => loadPosts(), [])
const { data: tagsData } = useLocalizedContent(() => loadTags(), [])

const siteAuthor = SITE.author
const currentTag = computed(() => (typeof route.query.tag === 'string' ? route.query.tag : ''))
const fromPath = computed(() => (typeof route.query.from === 'string' ? route.query.from : ''))
const filteredDocs = computed(() => {
  const tag = currentTag.value
  if (!tag) return postData.value
  return postData.value.filter((p) => p.tags.includes(tag))
})

interface TagItem {
  name: string
  count: number
  level: 'lg' | 'md' | 'sm'
}

const tagList = computed<TagItem[]>(() => {
  const map = new Map<string, number>()
  postData.value.forEach((p) => p.tags.forEach((t) => map.set(t, (map.get(t) || 0) + 1)))
  const list = Array.from(map.entries())
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => a.name.localeCompare(b.name))
  const total = list.length
  const third = Math.ceil(total / 3)
  const rank = new Map(
    [...list].sort((a, b) => b.count - a.count).map((item, idx) => [item.name, idx]),
  )
  return list.map((item) => {
    const idx = rank.get(item.name) ?? 0
    return { ...item, level: idx < third ? 'lg' : idx < third * 2 ? 'md' : 'sm' }
  })
})
const postCount = computed(() => postData.value.length)
const tagCount = computed(() => tagsData.value.length)
const totalWords = computed(() => postData.value.reduce((sum, p) => sum + p.wordCount, 0))
const totalWordsDisplay = computed(() => {
  const n = totalWords.value
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(n % 1_000_000 ? 1 : 0) + 'M'
  if (n >= 1_000) return (n / 1_000).toFixed(n % 1_000 ? 1 : 0) + 'K'
  return String(n)
})

function clearTag() {
  const q = { ...route.query }
  delete q.tag
  delete q.from
  q.page = '1'
  router.push({ path: route.path, query: q }).catch(() => {})
  scrollToTop()
}

function goBackFromTag() {
  if (fromPath.value) {
    router.push(fromPath.value).catch(() => {})
  }
}

function goTag(tag: string) {
  goToTag(router, tag, {
    locale: toSupportedLocale(locale.value),
    query: { ...route.query, page: '1' },
  })
}

watch(locale, (newLocale, oldLocale) => {
  if (newLocale !== oldLocale && currentTag.value) {
    clearTag()
  }
})
</script>

<style scoped>
.home-section {
  overflow: hidden;
}

.hero {
  position: relative;
  padding: var(--sp-24) 0 var(--sp-20);
  display: grid;
  grid-template-columns: minmax(0, 1fr) auto;
  align-items: end;
  gap: var(--sp-16);
}

.hero__main {
  min-width: 0;
}

.hero__greeting {
  display: inline-flex;
  align-items: center;
  gap: 10px;
  color: var(--fg-2);
  font-size: var(--text-sm);
  font-weight: 500;
  letter-spacing: 0.02em;
  margin-bottom: var(--sp-6);
}

.hero__greeting-mark {
  font-family: var(--font-mono);
  font-size: 13px;
  font-weight: 600;
  color: var(--primary);
}

.hero__name {
  font-size: clamp(46px, 7.5vw, 76px);
  font-weight: 700;
  letter-spacing: -0.03em;
  line-height: 1;
  margin-bottom: var(--sp-6);
  color: var(--fg);
}

.hero__bio {
  font-size: var(--text-lg);
  color: var(--fg-2);
  line-height: 1.65;
  max-width: 540px;
  margin: 0;
}

.hero__stats {
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: var(--sp-8);
  padding: 0 0 var(--sp-2);
  white-space: nowrap;
  min-width: 320px;
}

.hero__stat {
  display: flex;
  flex-direction: column;
  gap: 4px;
}

.hero__stat-num {
  font-size: 42px;
  font-weight: 700;
  letter-spacing: -0.02em;
  line-height: 1;
  color: var(--fg);
  font-variant-numeric: tabular-nums;
}

.hero__stat-label {
  font-size: 12px;
  color: var(--fg-3);
  letter-spacing: 0.08em;
  text-transform: uppercase;
  font-weight: 600;
}

.tags-row {
  display: flex;
  flex-wrap: wrap;
  gap: var(--sp-3);
  padding: var(--sp-8) 0 var(--sp-8);
}

.posts-header {
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  padding: var(--sp-8) 0 var(--sp-6);
  border-top: 1px solid var(--line);
  gap: var(--sp-3);
  flex-wrap: wrap;
}

.posts-header__title {
  font-size: var(--text-2xl);
  font-weight: 700;
  letter-spacing: -0.015em;
  color: var(--fg);
  margin: 0;
  display: flex;
  align-items: baseline;
  gap: var(--sp-3);
  flex-wrap: wrap;
}

.posts-header__title::before {
  content: '';
  width: 3px;
  height: 20px;
  border-radius: 99px;
  background: var(--primary);
  flex-shrink: 0;
  align-self: center;
}

.posts-header__title--tag::before {
  display: none;
}

.posts-header__actions {
  display: flex;
  align-items: center;
  gap: var(--sp-2);
  flex-wrap: wrap;
}

.posts-header__tag-title {
  color: var(--primary);
}

.chip-close {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  font-size: var(--text-sm);
  font-weight: 500;
  background: transparent;
  border: none;
  color: var(--fg-3);
  padding: 6px 12px;
  border-radius: var(--radius-pill);
  cursor: pointer;
  transition:
    color 0.14s ease,
    background-color 0.14s ease;
}

.chip-close:hover {
  color: var(--primary);
  background: var(--tint);
}

.chip-close i {
  font-size: 12px;
}

@media (max-width: 991px) {
  .hero {
    padding: var(--sp-16) 0 var(--sp-12);
    grid-template-columns: 1fr;
    gap: var(--sp-10);
  }

  .hero__name {
    font-size: clamp(42px, 9vw, 60px);
  }

  .hero__stats {
    justify-self: start;
  }
}

@media (max-width: 767px) {
  .hero {
    padding: var(--sp-12) 0 var(--sp-10);
  }

  .hero__name {
    font-size: 40px;
  }

  .hero__bio {
    font-size: 16px;
  }

  .hero__stats {
    gap: var(--sp-6);
  }

  .hero__stat-num {
    font-size: 32px;
  }
}
</style>
