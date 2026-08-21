<template>
  <div class="page-section category-view">
    <div class="category-head">
      <h1 class="article-title" v-reveal>{{ pageTitle }}</h1>
    </div>

    <div v-for="(category, index) in categoryList" :key="index" class="cat-section">
      <div class="cat-section__header" v-reveal>
        <h2 class="cat-section__title">
          <i
            :class="['fas', sectionIcon(category.title)]"
            class="cat-section__icon"
            aria-hidden="true"
          ></i>
          {{ category.title }}
        </h2>
        <span class="cat-section__count">{{ sectionCount(category) }}</span>
      </div>

      <div class="cat-grid">
        <div
          v-for="(item, idx) in category.items"
          :key="idx"
          class="cat-card"
          v-reveal
          :style="{ '--reveal-delay': Math.min(Number(idx), 5) * 40 + 'ms' }"
        >
          <div class="cat-card__head">
            <div class="cat-card__name">{{ item.title }}</div>
            <div class="cat-card__ext-links">
              <a
                v-if="item.github"
                :href="normalizeUrl(item.github)"
                target="_blank"
                rel="noopener noreferrer"
                class="cat-card__ext-link"
                :aria-label="'GitHub'"
                ><i class="fab fa-github"></i
              ></a>
              <a
                v-if="item.doi"
                :href="normalizeDoi(item.doi)"
                target="_blank"
                rel="noopener noreferrer"
                class="cat-card__ext-link"
                :aria-label="'DOI'"
                ><i class="fas fa-link"></i
              ></a>
            </div>
          </div>
          <p class="cat-card__desc">{{ item.desc }}</p>
          <div v-if="category.title === notesTitle" class="cat-card__stats">
            <span v-if="item.stats?.postsCount" class="cat-stat"
              ><i class="fas fa-file-lines"></i
              >{{ t('countPosts', { count: item.stats.postsCount }) }}</span
            >
            <span v-if="item.stats?.totalWords" class="cat-stat"
              ><i class="fas fa-font"></i
              >{{ t('countWords', { count: item.stats.totalWords }) }}</span
            >
            <span v-if="getLatestDate(item)" class="cat-stat"
              ><i class="fas fa-clock"></i>{{ getLatestDate(item) }}</span
            >
          </div>
          <div v-else class="cat-card__stats">
            <span v-if="item.language" class="cat-stat"
              ><i class="fas fa-code"></i>{{ item.language }}</span
            >
            <span v-if="item.year" class="cat-stat"
              ><i class="fas fa-calendar"></i>{{ item.year }}</span
            >
            <span v-if="item.license" class="cat-stat"
              ><i class="fas fa-scale-balanced"></i>{{ item.license }}</span
            >
          </div>

          <div v-if="Array.isArray(item.tags) && item.tags.length" class="cat-card__tags">
            <span v-for="(tag, tIdx) in item.tags" :key="tIdx" class="cat-card__tag">{{
              tag
            }}</span>
          </div>

          <div class="cat-card__links">
            <a
              v-if="hasExternalLink(item) || item.root"
              class="cat-card__link"
              @click.prevent="handleSeeMore(item)"
              >{{ seeMoreText }}<i class="fas fa-arrow-right"></i
            ></a>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup lang="ts">
defineOptions({ name: 'CategoryView' })
import { computed } from 'vue'
import { useI18n } from 'vue-i18n'
import { useHead } from '@unhead/vue'
import { useRouter } from 'vue-router'
import { useLocalizedContent } from '@/composables/useLocalizedContent'
import { toLocalePath } from '@/utils/navigation'
import { loadCategories } from '@/utils/contentLoader'
import type { CategoryItem, CategorySection } from '@/types'

useHead({ title: 'zorrooz’s blog - Categories' })

const { t } = useI18n()
const router = useRouter()

const { data: categoryList } = useLocalizedContent(() => loadCategories(), [])

const pageTitle = computed(() => t('categories'))
const notesTitle = computed(() => t('notes'))
const projectsTitle = computed(() => t('projects'))
const topicsTitle = computed(() => t('topics'))
const seeMoreText = computed(() => t('seeMore'))

function sectionIcon(title: string) {
  if (title === notesTitle.value) return 'fa-book-open'
  if (title === projectsTitle.value) return 'fa-folder-open'
  if (title === topicsTitle.value) return 'fa-flask'
  return 'fa-folder'
}

function sectionCount(category: CategorySection) {
  const items = category.items
  if (category.title === notesTitle.value) {
    const posts = items.reduce((s, it) => s + it.stats.postsCount, 0)
    return t('countPosts', { count: posts })
  }
  if (category.title === projectsTitle.value) return t('countProjects', { count: items.length })
  return t('countTopics', { count: items.length })
}

function getLatestDate(item: CategoryItem) {
  return item.stats.latestDate || ''
}

function normalizeUrl(value: string | undefined) {
  if (!value || !value.trim()) return ''
  if (/^https?:\/\//i.test(value)) return value
  return 'https://' + value.replace(/^\/+/, '')
}

function normalizeDoi(value: string | undefined) {
  if (!value || !value.trim()) return ''
  if (/^https?:\/\//i.test(value)) return value
  if (/^10\.\d{4,9}\//.test(value)) return 'https://doi.org/' + value
  return 'https://' + value.replace(/^\/+/, '')
}

function hasExternalLink(item: CategoryItem) {
  return !!(item.url || item.github || item.doi)
}

function handleSeeMore(item: CategoryItem) {
  if (item.root) {
    router.push(toLocalePath(item.root)).catch((err) => {
      if (err.name !== 'NavigationDuplicated' && !err.toString().includes('Navigation cancelled')) {
        console.error('Navigation error:', err)
      }
    })
    return
  }

  const fields: Array<'url' | 'github' | 'doi'> = ['url', 'github', 'doi']
  for (const field of fields) {
    const raw = item[field]
    if (!raw) continue
    const url = field === 'doi' ? normalizeDoi(raw) : normalizeUrl(raw)
    if (url) {
      window.open(url, '_blank', 'noopener,noreferrer')
      return
    }
  }
}
</script>

<style scoped>
.category-head {
  margin-bottom: var(--sp-16);
}

.cat-section {
  margin-bottom: var(--sp-24);
}

.cat-section__header {
  display: flex;
  align-items: center;
  gap: var(--sp-4);
  margin-bottom: var(--sp-8);
  padding-bottom: var(--sp-5);
  border-bottom: 1px solid var(--line);
}

.cat-section__title {
  font-size: var(--text-3xl);
  font-weight: 700;
  letter-spacing: -0.02em;
  color: var(--fg);
  display: flex;
  align-items: center;
  gap: var(--sp-3);
  margin: 0;
}

.cat-section__icon {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 30px;
  height: 30px;
  font-size: 13px;
  color: var(--primary);
  background: var(--tint);
  border-radius: var(--radius-btn);
  flex-shrink: 0;
}

/* 计数：文章页风格纯文本 meta（去胶囊底），数字等宽对齐 */
.cat-section__count {
  font-size: var(--text-xs);
  font-weight: 600;
  color: var(--fg-3);
  letter-spacing: 0.04em;
  font-variant-numeric: tabular-nums;
  margin-left: auto;
}

.cat-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(310px, 1fr));
  gap: var(--sp-5);
}

.cat-card {
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  padding: var(--sp-8);
  display: flex;
  flex-direction: column;
  position: relative;
  transition: border-color 0.18s ease;
}

.cat-card:hover {
  border-color: color-mix(in srgb, var(--primary) 30%, transparent);
}

.cat-card__head {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: var(--sp-3);
}

.cat-card__name {
  font-size: var(--text-lg);
  font-weight: 600;
  letter-spacing: -0.01em;
  color: var(--fg);
  line-height: 1.3;
  min-width: 0;
  transition: color 0.14s ease;
}

.cat-card:hover .cat-card__name {
  color: var(--primary);
}

.cat-card__ext-links {
  display: flex;
  align-items: center;
  gap: 4px;
  flex-shrink: 0;
}

.cat-card__ext-link {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 32px;
  height: 32px;
  border-radius: 50%;
  color: var(--fg-3);
  font-size: 14px;
  text-decoration: none;
  transition:
    color 0.14s ease,
    background-color 0.14s ease;
}

.cat-card__ext-link:hover {
  color: var(--primary);
  background: var(--tint);
}

.cat-card__desc {
  font-size: var(--text-md);
  color: var(--fg-2);
  line-height: 1.65;
  margin: var(--sp-4) 0 var(--sp-5);
}

.cat-card__stats {
  display: flex;
  gap: var(--sp-3);
  flex-wrap: wrap;
  align-items: center;
  margin-top: var(--sp-2);
}

.cat-stat {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  font-size: var(--text-xs);
  font-weight: 500;
  color: var(--fg-3);
  letter-spacing: 0.02em;
  font-variant-numeric: tabular-nums;
}

.cat-stat i {
  font-size: 10px;
  color: var(--primary);
}

.cat-card__tags {
  display: flex;
  flex-wrap: wrap;
  gap: var(--sp-2);
  margin-top: var(--sp-4);
}

.cat-card__tag {
  font-size: var(--text-xs);
  font-weight: 500;
  padding: 3px 11px;
  border-radius: var(--radius-pill);
  background: var(--surface-2);
  color: var(--fg-2);
}

.cat-card__links {
  display: flex;
  gap: var(--sp-5);
  margin-top: auto;
  padding-top: var(--sp-5);
  font-size: var(--text-sm);
}

.cat-card__link {
  color: var(--app-link);
  font-weight: 500;
  display: inline-flex;
  align-items: center;
  gap: 6px;
  cursor: pointer;
  background: none;
  border: none;
  padding: 0;
  font-size: var(--text-sm);
  transition:
    color 0.14s ease,
    gap 0.14s ease;
}

.cat-card__link i {
  font-size: 12px;
}

.cat-card__link:hover {
  color: var(--primary-strong, var(--primary));
  gap: 9px;
}

@media (max-width: 767px) {
  .cat-grid {
    grid-template-columns: 1fr;
  }

  .cat-card {
    padding: var(--sp-6);
  }
}
</style>
