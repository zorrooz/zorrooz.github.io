<template>
  <div class="page-section resource-view">
    <header class="resource-head">
      <h1 class="article-title">{{ pageTitle }}</h1>
      <p class="resource-subtitle">
        <i class="fas fa-circle-info resource-head__icon"></i>{{ pageSubtitle }}
      </p>
    </header>

    <div class="res-layout">
      <aside class="res-sidebar">
        <div v-for="category in resources" :key="category.title" class="res-group">
          <div class="res-group__label">
            <span>{{ category.title }}</span>
            <span class="res-group__count num">{{ category.children?.length || 0 }}</span>
          </div>
          <div class="res-group__items">
            <button
              v-for="sub in category.children"
              :key="sub.title"
              class="res-item"
              :class="{ 'res-item--active': isActiveSub(sub) }"
              @click="selectSub(sub)"
            >
              {{ sub.title }}
            </button>
          </div>
        </div>
      </aside>

      <div class="res-divider" aria-hidden="true"></div>

      <main class="res-main">
        <div class="res-groups">
          <section v-for="group in groups" :key="group.title" class="res-group-block">
            <h3 class="res-group-block__title">{{ group.title }}</h3>
            <div class="res-grid">
              <a
                v-for="item in group.items"
                :key="item.name"
                :href="item.url"
                target="_blank"
                rel="noopener noreferrer"
                class="res-card"
              >
                <div class="res-card__head">
                  <span class="res-card__name">{{ item.name }}</span>
                  <span class="res-card__ext-links">
                    <span v-if="isDoi(item.url)" class="res-card__ext-link" aria-label="DOI"
                      ><i class="fas fa-link"></i
                    ></span>
                  </span>
                </div>
                <p class="res-card__desc">{{ item.desc }}</p>
                <div class="res-card__footer">
                  <span class="res-card__url">{{ displayUrl(item.url) }}</span>
                  <i class="fas fa-arrow-up-right-from-square res-card__arrow"></i>
                </div>
              </a>
            </div>
          </section>
        </div>
      </main>
    </div>
  </div>
</template>

<script setup lang="ts">
defineOptions({ name: 'ResourceView' })
import { computed, ref, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useLocalizedContent } from '@/composables/useLocalizedContent'
import { loadResources } from '@/utils/contentLoader'
import type { ResourceNode } from '@/types'

const { t } = useI18n()

const { data: resources } = useLocalizedContent(() => loadResources(), [])
const activeSub = ref<ResourceNode | null>(null)

const pageTitle = computed(() => t('resources'))
const pageSubtitle = computed(() => t('resourceSubtitle'))

const groups = computed<ResourceNode[]>(() => {
  const sub = activeSub.value
  if (!sub) return []
  if (sub.children?.length) {
    return sub.children.filter((g) => g.items?.length)
  }
  return sub.items?.length ? [{ title: sub.title, items: sub.items }] : []
})

function selectSub(sub: ResourceNode) {
  activeSub.value = sub
}

function isActiveSub(sub: ResourceNode) {
  return activeSub.value === sub
}

function displayUrl(url: string) {
  if (!url) return ''
  return url.replace(/^https?:\/\//, '').replace(/\/$/, '')
}

function isDoi(url: string) {
  return !!url && (url.includes('doi.org') || /^10\.\d{4,9}\//.test(url))
}

watch(
  resources,
  (list) => {
    activeSub.value = list[0]?.children?.[0] || null
  },
  { immediate: true },
)
</script>

<style scoped>
.resource-head {
  margin-bottom: var(--sp-12);
}

.resource-head__icon {
  font-size: 14px;
  color: var(--primary);
  margin-right: 8px;
}

.resource-subtitle {
  color: var(--fg-2);
  font-size: var(--text-md);
  margin-top: var(--sp-3);
  margin-bottom: 0;
}

.res-layout {
  display: grid;
  grid-template-columns: 230px 1px 1fr;
  gap: 0;
  min-height: 420px;
}

.res-sidebar {
  padding: var(--sp-2) var(--sp-5) var(--sp-2) 0;
  position: sticky;
  top: 92px;
  align-self: start;
  max-height: calc(100vh - 116px);
  overflow-y: auto;
}

.res-group {
  margin-bottom: var(--sp-8);
}

.res-group__label {
  font-size: 12px;
  font-weight: 600;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--fg-3);
  padding: 0 var(--sp-2);
  margin-bottom: var(--sp-3);
  display: flex;
  align-items: center;
  gap: 8px;
}

.res-group__count {
  margin-left: auto;
  font-size: 11px;
  font-weight: 600;
  color: var(--fg-3);
  background: var(--surface-2);
  border-radius: 99px;
  padding: 1px 8px;
}

.res-group__items {
  display: flex;
  flex-direction: column;
}

.res-item {
  display: flex;
  align-items: center;
  gap: var(--sp-2);
  font-size: var(--text-base);
  color: var(--fg-2);
  padding: 8px 12px;
  border-radius: var(--radius-sm);
  cursor: pointer;
  transition:
    color var(--dur-fast) ease,
    background-color var(--dur-fast) ease;
  width: 100%;
  text-align: left;
  border: none;
  background: transparent;
  margin-bottom: 2px;
  font-weight: 500;
}

.res-item:hover {
  background: var(--surface-2);
  color: var(--fg);
}

/* 当前选中项 hover 保持主色（与 NavigationTree 约定一致），不落入 hover 的黑色 */
.res-item--active:hover {
  color: var(--primary);
  font-weight: 600;
  background: var(--tint-strong);
}

.res-item--active {
  color: var(--primary);
  font-weight: 600;
  background: var(--tint);
}

.res-divider {
  width: 1px;
  background: var(--line);
  align-self: stretch;
}

.res-main {
  padding: var(--sp-2) 0 var(--sp-2) var(--sp-10);
  min-width: 0;
}

.res-group-block {
  margin-bottom: var(--sp-8);
}

.res-group-block:last-child {
  margin-bottom: 0;
}

.res-group-block__title {
  font-size: var(--text-xl);
  font-weight: 700;
  letter-spacing: -0.015em;
  color: var(--fg);
  margin: 0 0 var(--sp-4);
  padding-bottom: var(--sp-3);
  border-bottom: 1px solid var(--line);
  display: flex;
  align-items: center;
  gap: var(--sp-3);
}

.res-group-block__title::before {
  content: '';
  width: 3px;
  height: 18px;
  border-radius: 99px;
  background: var(--primary);
  flex-shrink: 0;
}

.res-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
  gap: var(--sp-4);
}

.res-card {
  display: flex;
  flex-direction: column;
  padding: var(--sp-6);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  background: var(--surface);
  text-decoration: none;
  color: var(--fg);
  transition:
    border-color 0.18s ease,
    box-shadow 0.18s ease;
}

.res-card:hover {
  border-color: color-mix(in srgb, var(--primary) 30%, transparent);
  box-shadow: var(--shadow-soft);
}

.res-card__head {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: var(--sp-3);
}

.res-card__name {
  font-size: var(--text-lg);
  font-weight: 600;
  color: var(--fg);
  min-width: 0;
  transition: color 0.14s ease;
}

.res-card:hover .res-card__name {
  color: var(--primary);
}

.res-card__ext-links {
  display: flex;
  align-items: center;
  gap: 4px;
  flex-shrink: 0;
}

.res-card__ext-link {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 30px;
  height: 30px;
  border-radius: 50%;
  color: var(--primary-muted);
  font-size: 14px;
  transition:
    color 0.14s ease,
    background-color 0.14s ease;
}

.res-card__ext-link:hover {
  color: var(--primary);
  background: var(--tint);
}

.res-card__ext-link i {
  pointer-events: none;
}

.res-card__desc {
  font-size: var(--text-base);
  color: var(--fg-2);
  line-height: 1.6;
  margin: var(--sp-3) 0 var(--sp-4);
  flex-grow: 1;
}

.res-card__footer {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: var(--sp-3);
  border-top: 1px solid var(--line);
  padding-top: var(--sp-3);
}

.res-card__url {
  font-size: var(--text-xs);
  color: var(--fg-3);
  font-family: var(--font-mono);
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.res-card__arrow {
  font-size: 12px;
  color: var(--fg-3);
  padding: 4px;
  opacity: 0;
  transform: translateX(-4px);
  transition:
    opacity var(--dur-base) ease,
    transform var(--dur-base) ease,
    color var(--dur-fast) ease;
}

.res-card:hover .res-card__arrow {
  opacity: 1;
  transform: translateX(0);
}

.res-card__arrow:hover {
  color: var(--primary);
}

@media (max-width: 1023px) {
  .res-layout {
    grid-template-columns: 1fr;
  }

  .res-divider {
    display: none;
  }

  .res-sidebar {
    position: static;
    max-height: none;
    overflow-y: visible;
    padding: 0;
    display: flex;
    flex-direction: column;
    gap: var(--sp-5);
    margin-bottom: var(--sp-8);
  }

  .res-group {
    margin-bottom: 0;
  }

  .res-group__items {
    flex-direction: row;
    gap: var(--sp-2);
    overflow-x: auto;
    padding-bottom: 4px;
    -webkit-overflow-scrolling: touch;
    scrollbar-width: none;
    -ms-overflow-style: none;
  }

  .res-group__items::-webkit-scrollbar {
    display: none;
  }

  .res-item {
    flex-shrink: 0;
    width: auto;
    margin-bottom: 0;
    border-radius: var(--radius-pill);
    background: var(--surface-2);
    padding: 7px 15px;
    font-size: var(--text-sm);
  }

  .res-item--active {
    background: var(--tint);
  }

  .res-main {
    padding-left: 0;
  }
}
</style>
